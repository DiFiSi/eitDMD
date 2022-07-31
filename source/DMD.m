classdef DMD < handle
    %DMD Summary of this class goes here
    %   Detailed explanation goes here
    
    % TODO
%     Try running DMD at each sample
%     How long is a prediction available for?
%     Define lengths of windows
    
    properties(Access = private)
       PerfCm = buildcmap('yrk');
       BarCm = buildcmap('bwr');
    end
    
    properties(Constant)
       ProjErrThresh = 0.5;
       Rmin = 30;
       Rmax = 150;
    end
    
    properties(Access = public)
       Rank
        
       Ux
       Sx
       Vx
       
       Uy
       Sy
       Vy
      
       Pinvx
       Pinvy
       Q
       
       Atilde
       lambda
       omega
       Phi
       b
       Norms
       Mags
       Freqs
    end
    
    methods
        function obj = DMD()
            %DMD Construct an instance of this class
        end
        
        % Compressing data for compressed sensing DMD (according to DMD Book)
        function [cX, cY, C] = compressData(obj, X, Y, m, projType)
            data = [X(:,1),Y];
            n = size(data, 1);
            
            % Create random measurement matrix
            C = zeros(m, n);
%             Theta = zeros(m, n);
            xMeas = zeros(n, 1);
            for i = 1:m
                xMeas = 0 * xMeas;
                if(projType == 1)
                    xMeas = rand(n,1);
                    
                elseif(projType == 2)
                    xMeas = randn(n,1);
                    
                elseif(projType == 3)
                    xMeas(ceil(n * rand),ceil(n * rand)) = 1;
                end
                
                C(i,:) = reshape(xMeas, n, 1);
%                 Theta(i,:) = reshape((ifft2(xMeas)), 1, n);
            end
            
            data = C * data;
            
            cX = data(:,1:end -1);
            cY = data(:,2:end);
        end
        
        function Xout = hankelTransform(obj, X, maxDelay)
            if maxDelay == 0
                Xout = X;
            else
                [nStates,nSnaps] = size(X);
                Xout = zeros(nStates * maxDelay, nSnaps - maxDelay);
                for i = 0:maxDelay - 1
                    startRow = nStates * i + 1;
                    finishRow = nStates * (i + 1);

                    startCol = 1 + i;
                    finishCol = nSnaps - (maxDelay - i);

                    Xout(startRow:finishRow,:) = X(:,startCol:finishCol);
                end
            end
        end
        
        % Getting Atilde
        function fitKoopman(obj, X, Y, rank, tSize, tp) % TODO: Faster truncated SVD, Optimal SVD truncation for EIT
            if ~exist('tp','var')
               tp = 'standard'; 
            end
            
            [obj.Ux,obj.Sx,obj.Vx,obj.Rank] = truncSVD(X,rank,tSize);
            
            switch(tp)
                case 'stream'
                    Xtilde = obj.Ux' * X;
                    
                    [obj.Uy,obj.Sy,obj.Vy] = truncSVD(Y,obj.Rank,tSize);
                    Ytilde = obj.Uy' * Y;
                    
                    obj.Q = Ytilde * Xtilde';
                    obj.Pinvx = Xtilde * Xtilde';
                    obj.Pinvy = Ytilde * Ytilde';
                    
                    obj.Atilde = obj.Ux' * obj.Uy * obj.Q / obj.Pinvx;
                    
                case 'fb'
                    forwAtilde = obj.Ux' * Y * obj.Vx / obj.Sx;

                    [obj.Uy,obj.Sy,obj.Vy] = truncSVD(Y,obj.Rank,tSize);
                    backAtilde = obj.Uy' * X * obj.Vy / obj.Sy;

                    obj.Atilde = (forwAtilde / backAtilde) ^ 0.5;
                    
                case 'tls'
                    Z = [X;Y];
                    
                    [Uz, ~, ~] = svd(Z, 'econ');

                    U11 = Uz(1:obj.Rank, 1:obj.Rank);
                    U21 = Uz(obj.Rank + 1:end, 1:obj.Rank);
                    
                    obj.Atilde = U21 / U11;
                    
                case 'standard'
                    obj.Atilde = obj.Ux' * Y * obj.Vx / obj.Sx;
            end
        end
        
        function updateStream(obj, x, y, rho)
            % Get the 2-norm of the original vectors
            normx = norm(x,2);
            normy = norm(y,2);
            
            % Project new vectors onto singular vectors of data
            xtilde = obj.Ux' * x;
            ytilde = obj.Uy' * y;
            
            % Calculate error of projection
            ex = x - obj.Ux * xtilde;
            normex = norm(ex,2);
            ey = y - obj.Uy * ytilde;
            normey = norm(ey,2);
            
            % Increase rank of Ux and Uy if projection error high
            if normex / normx > obj.ProjErrThresh
                % Get new vector for Ux matrix from normalized error
                ux = ex / normex;
                obj.Ux = [obj.Ux,ux];
                
                % Expand Pinv and Q
                obj.Pinvx = [[obj.Pinvx,zeros(size(obj.Pinvx,1),1)];zeros(1,size(obj.Pinvx,2))];
                obj.Q = [obj.Q,zeros(size(obj.Q,1),1)];
            end
            
            if normey / normy > obj.ProjErrThresh
                % Get new vector for Uy matrix from normalized error
                uy = ey / normey;
                obj.Uy = [obj.Uy,uy];
                
                % Expand Pinv and Q
                obj.Pinvy = [[obj.Pinvy,zeros(size(obj.Pinvy,1),1)];zeros(1,size(obj.Pinvy,2))];
                obj.Q = [obj.Q,zeros(size(obj.Q,1),1)];
            end
            
            % Reduce rank of Ux and Uy if current rank high
            if size(obj.Ux,1) > obj.Rmax
               [eigVecs,eigVals] = eig(obj.Pinvx,'vector'); 
               [eigVals,idx] = sort(eigVals,'descend');
               qx = eigVecs(:,idx(1:obj.Rmin));
               
               obj.Ux = obj.Ux * qx;
               obj.Q = obj.Q * qx;
               obj.Pinvx = diag(eigVals(1:obj.Rmin));
            end
            
            if size(obj.Uy,1) > obj.Rmax
               [eigVecs,eigVals] = eig(obj.Pinvy,'vector'); 
               [eigVals,idx] = sort(eigVals,'descend');
               qy = eigVecs(:,idx(1:obj.Rmin));
               
               obj.Uy = obj.Uy * qy;
               obj.Q = qy' * obj.Q;
               obj.Pinvy = diag(eigVals(1:obj.Rmin));
            end
            
            xtilde = obj.Ux' * x;
            ytilde = obj.Uy' * y;
            
            % Update regression
            obj.Q = rho * obj.Q + ytilde * xtilde';
            obj.Pinvx = rho * obj.Pinvx + xtilde * xtilde';
            obj.Pinvy = rho * obj.Pinvy + ytilde * ytilde';
            obj.Atilde = obj.Ux' * obj.Uy * obj.Q / obj.Pinvx;
        end
        
        % Find Phi, b, and lambda
        function fitModes(obj, X, Y, fs, tp, fX, fY)
            % Get eigenvalues and eigenvectors
            [W,D] = eig(obj.Atilde);
            obj.lambda = diag(D);
            obj.omega = log(obj.lambda)/(1/fs);
            
            % Only used in case of compressed sensing
            if nargin > 9
               Y = fY;
               X = fX;
            end
            
            fSize = size(X, 1);
            obj.Rank = length(obj.lambda);
            
            % Calculate DMD modes and amplitudes (interpreted in absolute)
            switch(tp)
                case 'exact-mod'
                    tmp = Y * obj.Vx / obj.Sx;
                    obj.Phi = zeros(fSize, obj.Rank);
                    for i = 1:obj.Rank
                        obj.Phi(:,i) = 1/obj.lambda(i) * tmp * W(:,i);
                    end
                    obj.b = inv(D) / obj.Phi * X(:,2);
                    
                % Exact DMD from Peter J. Schmid (https://www.youtube.com/watch?v=xAYimi7x4Lc)
                case 'exact'
                    obj.Phi = Y * obj.Vx / obj.Sx * W;
                    obj.b = obj.Phi \ X;
                    
                % Simple vanilla DMD from DMD Book
                case 'vanilla'
                    obj.Phi = obj.Ux * W;
                    obj.b = obj.Phi \ X;
                    
                case 'hankel'
                    obj.Phi = Y * obj.Vx / obj.Sx * W;
                    obj.b = obj.Phi \ X(:,1);
                    
                otherwise
                    error('Invalid fitting type.');
            end
        end
        
        % Reconstruct data from modes
        function [Xcomps, Xfinal] = reconData(obj, X, tp)
            if ~exist('tp','var')
               tp = 'standard';
            end
%             rank = length(obj.lambda);
            fSize = size(obj.Phi, 1);
%             t = (0:tSize - 1) * (1 / fs);
            
            nPred = 1:size(X,2);
            nSnaps = length(nPred);
            Xcomps = cell(obj.Rank,1);
            Xfinal = zeros(fSize, nSnaps);
            
            switch(tp)
                case 'standard'
                    % Reconstruct data (discrete)
                    for k = 1:nSnaps
                        n = nPred(k);
                        for j = 1:obj.Rank
                            tmp = obj.lambda(j)^(n-1) * obj.b(j) * obj.Phi(:,j);
                            Xcomps{j}(:,k) = tmp;
                            Xfinal(:,k) = Xfinal(:,k) + tmp;
                        end
                    end
                    
                case 'stream'
                    Xfinal = obj.Uy * obj.Atilde * obj.Ux' * X;
                    Xcomps = Xfinal;
            end
            
%             switch(tp)
%                 % Discrete version
%                 
%                 case 'disc'
%                     % Reconstruct data (discrete)
%                     for k = 1:tSize
%                         for j = 1:obj.Rank
%                             tmp = obj.lambda(j)^(k-1) * obj.b(j) * obj.Phi(:,j);
%                             Xcomps{j}(:,k) = tmp;
%                             Xfinal(:,k) = Xfinal(:,k) + tmp;
%                         end
%                     end
%                     
%                 % Continuous version
%                 case 'cont'
%                     % Reconstruct data (continuous)
%                     for k = 1:tSize
%                         for j = 1:obj.Rank
%                             tmp = obj.Phi(:,j) * obj.b(j) * exp(obj.omega(j) * t(k));
%                             Xcomps{j}(:,k) = tmp;
%                             Xfinal(:,k) = Xfinal(:,k) + tmp;
%                         end
%                     end
%                     
%                 otherwise
%                     error('Invalid reconstruction type.');
%             end
        end
        
%         function error = reconError(obj, Xfinal, X)
%             pass;
%         end
        
        % Remove complex duplicates from modes
        function cleanModes(obj)
            uniqueIdx = (imag(obj.lambda) >= 0);
            obj.Rank = sum(uniqueIdx);
            obj.lambda = obj.lambda(uniqueIdx);
            obj.omega = obj.omega(uniqueIdx);
            obj.Phi = obj.Phi(:, uniqueIdx);
            obj.b = obj.b(uniqueIdx);
        end
        
        % Interpret modes into their norms, magnitudes and frequencies
        % (according to https://doi.org/10.1016/j.visinf.2021.06.003)
        function interpModes(obj)
            obj.Norms = abs(obj.b'.* vecnorm(obj.Phi,2,1));
            obj.Mags = abs(obj.lambda);
            obj.Freqs = abs(imag(obj.omega/2/pi));
        end
        
        % Plot half polar plot
        function damFreq(obj)
            % magnitude vs frequency plot
            obj.Norms = obj.Norms / max(obj.Norms);
            freqRef = 0:(1.1 * max(obj.Freqs));
            
            figure; 
            
            plot(freqRef, ones(size(freqRef)), 'LineStyle','--','Color','k'); 
            hold on; 
            scatter(obj.Freqs, obj.Mags, 50 * ones(size(obj.Freqs)), obj.Norms, 'filled', 'Marker','o');
            
            ylabel('Mode Damping Ratio [1]', 'Interpreter', 'Latex', 'FontSize', 12);
            xlabel('Mode Frequency [Hz]', 'Interpreter', 'Latex', 'FontSize', 12);
            cb = colorbar; colormap(obj.PerfCm); 
            ylabel(cb,'Mode Intensity [n.u.]', 'Interpreter', 'Latex', 'FontSize', 12);
            ylim([0 1.1 * max(obj.Mags)]); xlim([0 1.1 * max(obj.Freqs)]);
        end
        
        % Plot dominance structure (according to https://doi.org/10.1016/j.visinf.2021.06.003)
        function domStruct(obj, tSize)
%             rank = length(obj.lambda);
            
            % Dominance structure plot
            nDivs = 8;
            cm = obj.BarCm(round(linspace(1,256,nDivs - 1)), :);

            dom = zeros(obj.Rank, tSize);
            for k = 1:tSize
                for j = 1:obj.Rank
                    dom(j,k) = norm(obj.lambda(j) ^ (k - 1) * obj.b(j) * obj.Phi(:,j));
                end
            end
            
            divs = round(linspace(1, tSize, nDivs));
            finalDom = zeros(obj.Rank, nDivs);
            for j = 1:obj.Rank
                for i = 1:nDivs - 1
                    finalDom(j, i) = sum(dom(j, divs(i):divs(i + 1)), 2) ./ (divs(i + 1) - divs(i));
                end
            end

            [~, idx] = sort(obj.Freqs);
            
            figure; 
            
            ba = bar(finalDom(idx,:),'stacked', 'FaceColor','flat');
            for i = 1:nDivs - 1
                ba(i).CData = cm(i,:);
                ba(i).BarWidth = 1;
            end
            ylabel('Mode Intensity [n.u.]', 'Interpreter', 'Latex', 'FontSize', 12);
            xlabel('Mode Frequency [Ascending Ordinals]', 'Interpreter', 'Latex', 'FontSize', 12);
            
            colormap(cm); cb = colorbar; 
            ylabel(cb,'bins [1]', 'Interpreter', 'Latex', 'FontSize', 12);
            caxis([divs(1), divs(end)]);
            set(cb, 'Ticks', divs(1:end),'TickLabelInterpreter','latex')
        end
        
        
        function nyqPlot(obj)
            % Nyquist Plot
            norms = obj.Norms / max(obj.Norms);
            freqRef = 0:0.1:pi;

            figure;
            
            polarscatter(phase(obj.lambda), abs(obj.lambda), 50 * ones(size(obj.Freqs)), norms, 'filled');
            
            thetalim([0 180]); rlim([0 1.1 * max(obj.Mags)]);
            cb = colorbar; colormap(obj.PerfCm); 
            ylabel(cb,'Mode Intensity [n.u.]', 'Interpreter', 'Latex', 'FontSize', 12);
            title('Half Complex Plane', 'Interpreter', 'Latex', 'FontSize', 12);
            hold on; polarplot(freqRef, ones(size(freqRef)), 'k--');
        end
        
        % Power spectrum of signals against of modes
        function powSpect(obj, X, fs)
            [amp, ~, f] = ampPhaseFFT(X, fs);
            
            figure; 
            
            yyaxis left;
            maxAmps = repmat(max(amp,[],1,'omitnan'),size(amp,1),1);
            dbAmp = 20 * log10(amp./maxAmps);
            plot(f,mean(dbAmp,2,'omitnan'),'LineWidth',2,'Color','r');
            xlabel('Frequency [Hz]','Interpreter', 'Latex', 'FontSize', 12);
            ylabel('Amplitude [dB]', 'FontName', 'Arial', 'FontSize', 9);
            axis tight;
            grid on;
            xlim([0,10]);

            yyaxis right;
            stem(obj.Freqs,db(obj.Norms),'LineWidth', 1.5, 'LineStyle', '-',...
                     'Color', 'black',...
                     'MarkerFaceColor','black');
            hold on;
            scatter(obj.Freqs, db(obj.Norms), 50 * ones(size(obj.Freqs)), obj.Mags, 'filled', 'Marker','o',...
                'MarkerEdgeColor','black', "LineWidth",1);
            ylabel('Magnitude [dB]', 'FontName', 'Arial', 'FontSize', 9);
            xlabel('Frequency [Hz]', 'FontName', 'Arial', 'FontSize', 9);
            grid on;
            xlim([0,10]);
        end
        
        % Multi-resolution DMD (according to DMD book)
        function [lambdaCum, omegaCum, PhiCum, bCum] = mrFit(obj, X, Y, rank, nLevels, tSize, fX, fY)
%             nModes = sum(2 .^ [0:nLevels - 1]) * rank;
%             fSize = size(X, 1);
            
            lambdaCum = [];
            omegaCum = [];
            PhiCum = [];
            bCum = [];
            
            for l = 1:nLevels
                nBounds = 2 ^ (l - 1);
                winLen = floor(tSize / nBounds);
                
                finish = 0;
                for w = 0:nBounds - 1
                    start = finish + 1;
                    finish = start + winLen - 1;
                    
                    x = X(:, start:finish);
                    y = Y(:, start:finish);
                    
                    if nargin > 6
                        fx = fX(:, start:finish);
                        fy = fY(:, start:finish);
                    end
                    
                    if rank > winLen 
                        rank = winLen;
                    end
                    
                    fitKoopman(obj, x, y, rank, tSize);
                    
                    if nargin > 6
                        fitModes(obj, obj.Ur, obj.Sr, obj.Vr, obj.Atilde, x, y, 50, 'exact-mod', fx, fy);
                    else
                        fitModes(obj, obj.Ur, obj.Sr, obj.Vr,obj. Atilde, x, y, 50, 'exact-mod');
                    end
                    
                    obj.lambda = [obj.lambda; obj.lambda];
                    obj.omega = [obj.omega; obj.omega];
                    obj.Phi = [obj.Phi, obj.Phi];
                    obj.b = [obj.b; obj.b];
                end
            end
            
            cleanModes(obj);
        end
        
        % TODO: mode selection for EIT and improve mode selection for b
        % Selection of modes based on frequency, norm and magnitude ranges
        function boxSelection(obj, limsNorms, limsMags, limsFreqs)
            nnorms = obj.Norms / max(obj.Norms);
            
            isNorms = nnorms >= limsNorms(1) & nnorms <= limsNorms(2);
            isMags = obj.Mags >= limsMags(1) & obj.Mags <= limsMags(2);
            isFreqs = obj.Freqs >= limsFreqs(1) & obj.Freqs <= limsFreqs(2);
            
            isSelect = isNorms' & isMags & isFreqs;
            
            obj.Norms = obj.Norms(isSelect);
            obj.Mags = obj.Mags(isSelect);
            obj.Freqs = obj.Freqs(isSelect);
            
            obj.lambda = obj.lambda(isSelect);
            obj.omega = obj.omega(isSelect);
            obj.Phi = obj.Phi(:, isSelect);
            obj.b = obj.b(isSelect);
        end
        
        % Selection of modes based on ROI similarity between reconstruction
        % and singular vectors of data
        function [norms, mags, freqs, lambda, omega, Phi, b] = roiSelection(obj, snrThresh, norms, mags, freqs, lambda, omega, Phi, b, noiseMask, heartMask, lungMask)
            nModes = length(freqs);
            
            energy = zeros(nModes, 3);
            for m = 1:nModes
                phi = abs(reshape(real(Phi(:, m)), 32, 32));
                totalEnergy = sum(phi(:));
                
                tmp = phi .* heartMask;
                energy(m,1) = sum(tmp(:))/totalEnergy;
                
                tmp = phi .* lungMask;
                energy(m,2) = sum(tmp(:))/totalEnergy;
                
                tmp = phi .* noiseMask;
                energy(m,3) = sum(tmp(:))/totalEnergy;
            end 
            
            isSelect = (sum(energy(:,1:2), 2) ./ sum(energy, 2)) > snrThresh;
            
            norms = norms(isSelect);
            mags = mags(isSelect);
            freqs = freqs(isSelect);
            
            lambda = lambda(isSelect);
            omega = omega(isSelect);
            Phi = Phi(:, isSelect);
            b = b(isSelect); 
        end
        
        % TODO: reconstruction error
        
        % TODO: harmonical clustering implementation (after box selection)
        
        % TODO: sparsity promoting
        
        % UNDER CONSTRUCTION
%         function P = getP(obj, X)
%             P = pinv(X * X');
% %             [U,S,V] = svd(X * X','econ');
% %             P = V * diag(1/diag(S)) * U';
%         end
%         
%         % Trying to update Koopman directly (best to update SVD instead)
%         function [Anext, Pnext] = updateKoopman(obj, A, P, Xnext, Ynext, w)
%             gamma = 1 / (1 + Xnext' * P * Xnext);
%             Pnext = (1 / w) * P - gamma * P * (Xnext * Xnext') * P;
%             Anext = A + gamma * (Ynext - A * Xnext) * Xnext' * P;
%         end
%         
%         % Trying the weighted windowed online DMD from arXiv:1908.01047v3
%         % PROBLEM: THIS VERSION DOESN'T ALLOW CACLULATION OF B, WHICH MEANS
%         % THE MODES CAN#T BE RECONSTRUCTED AFTERWARDS - NO GOOD
%         function [Anext, Unext, Snext, Vnext] = updateKoopmanSVD(obj, A, X, Y, Xnext, Ynext, U, S, V)
%             % We start with UXk, SXk, VXk
%             % X큝 = Xnew(:,1:end-1);
%             % We can then calculate UX큝, SX큝, VX큝
%             % Xk+1 = Xnew;
%             z = [1, zeros(1, size(V,1) - 1)];
%             
%             % Step window forward
%             Xnew = [X(:,2:end), Xnext];
%             Ynew = [Y(:,2:end), Ynext];
%             
%             temp1 = S - U' * X(:,1) * z * V';
%             [Utemp1,Stemp1,Vtemp1] = svd(temp1,'econ'); % Utemp = Us큝; Stemp = Ss큝; Vtemp = Vs큝; 
%             
%             Upres = U * Utemp1; % Upres = UX큝;
%             Spres = Stemp1; % Spres = SX큝;
%             Vpres = V(:,2:end)' * Vtemp1'; % Vpres = VX큝; %wrong
%             
%             temp2 = [Spres, Upres' * Xnew];
%             [Utemp2, Stemp2, Vtemp2] = svd(temp2, 'econ');
%             
%             Unext = Upres * Utemp2;
%             Snext = Stemp2;
%             Vnext = [Vpres' * Vtemp2(:,1:end-1)'; Vtemp2(:,end)']; % wrong
%             
%             Anext = A + (Ynext - A * Xnext) * Vtemp2(end,:) / inv(Snext) * Unext'; % can I use Xnext = Unext * Snext * Vnext'?
%         end
    end
end

