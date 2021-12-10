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
    
    methods
        function obj = DMD()
            %DMD Construct an instance of this class
        end
        
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
        
        function [Ur, Sr, Vr, Atilde] = fitKoopman(obj, X, Y, rank, tSize) % TODO: Faster truncated SVD, Optimal SVD truncation for EIT
            fSize = size(X, 1);
            
            [U,S,V] = svd(X,'econ');

            if isempty(rank) || rank <= 0
                [U,S,V] = svd(X,'econ');
                beta = tSize / fSize;
                sigma = diag(S);
                tau = optimal_SVHT_coef(beta,0) * median(sigma);
                
                rank = find(sigma > tau,1,'last');
            end
            
            Ur = U(:,1:rank);
            Sr = S(1:rank,1:rank);
            Vr = V(:,1:rank);
            
            % Approximate low-rank Koopman operator
            Atilde = Ur' * Y * Vr / Sr;
        end
        
        function [lambda, omega, Phi, b] = fitModes(obj, Ur, Sr, Vr, Atilde, X, Y, fs, tp, fX, fY)
            % Get eigenvalues and eigenvectors
            [W,D] = eig(Atilde);
            lambda = diag(D);
            omega = log(lambda)/(1/fs);
            
            if nargin > 9
               Y = fY;
               X = fX;
            end
            
            fSize = size(X, 1);
            rank = length(lambda);
            
            % Calculate DMD modes and amplitudes (seen in absolute)
            switch(tp)
                case 'exact-mod'
                    tmp = Y * Vr / Sr;
                    Phi = zeros(fSize, rank);
                    for i = 1:rank
                        Phi(:,i) = 1/lambda(i) * tmp * W(:,i);
                    end
                    b = inv(D) / Phi * X(:,2);
                    
                case 'exact'
                    Phi = X' * Vr / Sr * W;
                    b = inv(Phi) \ X;
                    
                case 'vanilla'
                    Phi = Ur * W;
                    b = inv(Phi) \ X;
                    
                otherwise
                    error('Invalid fitting type.');
            end
        end
        
        function [Xcomps, Xfinal] = reconData(obj, lambda, omega, Phi, b, tSize, fs, tp)
            rank = length(lambda);
            fSize = size(Phi, 1);
            t = (0:tSize - 1) * (1 / fs);
            
            Xcomps = cell(rank,1);
            Xfinal = zeros(fSize, tSize);
            
            switch(tp)
                case 'disc'
                    % Reconstruct data (discrete)
                    for k = 1:tSize
                        for j = 1:rank
                            tmp = lambda(j)^(k-1) * b(j) * Phi(:,j);
                            Xcomps{j}(:,k) = tmp;
                            Xfinal(:,k) = Xfinal(:,k) + tmp;
                        end
                    end
                    
                case 'cont'
                    % Reconstruct data (continuous)
                    for k = 1:tSize
                        for j = 1:rank
                            tmp =  Phi(:,j) * b(j) * exp(omega(j) * t(k));
                            Xcomps{j}(:,k) = tmp;
                            Xfinal(:,k) = Xfinal(:,k) + tmp;
                        end
                    end
                    
                otherwise
                    error('Invalid reconstruction type.');
            end
        end
        
        function [lambda, omega, Phi, b, rank] = cleanModes(obj, lambda, omega, Phi, b)
            uniqueIdx = (imag(lambda) >= 0);
            rank = sum(uniqueIdx);
            lambda = lambda(uniqueIdx);
            omega = omega(uniqueIdx);
            Phi = Phi(:, uniqueIdx);
            b = b(uniqueIdx);
        end
        
        function [norms, mags, freqs] = interpModes(obj, lambda, Phi, b)
            norms = abs(b' .* vecnorm(Phi,2,1));
            mags = abs(lambda);
            freqs = 180 * angle(lambda)/(2 * pi ^ 2);
        end
        
        function magFreq(obj, norms, mags, freqs)
            % magnitude vs frequency plot
            norms = norms / max(norms);
            freqRef = 0:(1.1 * max(freqs));
            figure; 
            plot(freqRef, ones(size(freqRef)), 'LineStyle','--','Color','k'); 
            hold on; 
            scatter(freqs, mags, 50 * ones(size(freqs)), norms, 'filled', 'Marker','o');
            ylabel('Mode magnitude [1]', 'Interpreter', 'Latex', 'FontSize', 12);
            xlabel('Mode Frequency [Ascending Ordinals]', 'Interpreter', 'Latex', 'FontSize', 12);
            cb = colorbar; colormap(obj.PerfCm); 
            ylabel(cb,'Mode Norm [n.u.]', 'Interpreter', 'Latex', 'FontSize', 12);
            ylim([0 1.1 * max(mags)]); xlim([0 1.1 * max(freqs)]);
        end
        
        function domStruct(obj, freqs, lambda, Phi, b, tSize)
            rank = length(lambda);
            
            % Dominance structure plot
            nDivs = 8;
            cm = obj.BarCm(round(linspace(1,256,nDivs)), :);

            dom = zeros(rank, tSize);
            for k = 1:tSize
                for j = 1:rank
                    dom(j,k) = norm(lambda(j) ^ (k - 1) * b(j) * Phi(:,j));
                end
            end
            
            divs = round(linspace(1, tSize, nDivs));
            finalDom = zeros(rank, nDivs);
            for j = 1:rank
                for i = 1:nDivs - 1
                    finalDom(j, i) = sum(dom(j, divs(i):divs(i + 1)), 2) ./ (divs(i + 1) - divs(i));
                end
            end

            [~, idx] = sort(freqs);
            figure; 
            ba = bar(finalDom(idx,:),'stacked', 'FaceColor','flat');
            for i = 1:nDivs - 1
                ba(i).CData = cm(i,:);
                ba(i).BarWidth = 1;
            end
            ylabel('Mode Norm [n.u.]', 'Interpreter', 'Latex', 'FontSize', 12);
            xlabel('Mode Frequency [Ascending Ordinals]', 'Interpreter', 'Latex', 'FontSize', 12);
            
            colormap(cm); cb = colorbar; 
            ylabel(cb,'bins [1]', 'Interpreter', 'Latex', 'FontSize', 12);
            caxis([divs(1), divs(end)]);
            set(cb, 'Ticks', divs(1:end),'TickLabelInterpreter','latex')
        end
        
        function nyqPlot(obj, norms, freqs, mags, lambda)
            % Nyquist Plot
            norms = norms / max(norms);
            freqRef = 0:0.1:pi;

            figure; 
            polarscatter(phase(lambda), abs(lambda), 50 * ones(size(freqs)), norms, 'filled');
            thetalim([0 180]); rlim([0 1.1 * max(mags)]);
            cb = colorbar; colormap(obj.PerfCm); 
            ylabel(cb,'Mode Norm [n.u.]', 'Interpreter', 'Latex', 'FontSize', 12);
            title('Half Complex Plane', 'Interpreter', 'Latex', 'FontSize', 12);
            hold on; polarplot(freqRef, ones(size(freqRef)), 'k--');
        end
        
        function powSpect(obj, X, freqs, norms, fs)
            [pxx,f] = pwelch(X, [], [], [], fs);
            
            figure; 
            yyaxis left;
            shadedErrorBar(f,pow2db(pxx'),{@mean,@std},'lineProps',{'markerfacecolor','k'});
            xlabel('Frequency [Hz]','Interpreter', 'Latex', 'FontSize', 12);
            ylabel('PSD [dB/Hz]','Interpreter', 'Latex', 'FontSize', 12);
            axis tight;
            
            yyaxis right; 
            stem(freqs,norms,'LineWidth', 1.5, 'LineStyle', ':',...
                     'Color', 'red',...
                     'MarkerFaceColor','black',...
                     'MarkerEdgeColor','red', 'MarkerSize', 8);
            ylabel('Mode Norm [n.u.]','Interpreter', 'Latex', 'FontSize', 12);
        end
        
        function [lambdaCum, omegaCum, PhiCum, bCum] = mrFit(obj, X, Y, rank, nLevels, tSize)
%             nModes = sum(2 .^ [0:nLevels - 1]) * rank;
%             fSize = size(X, 1);
            
            [cX, cY, ~] = compressData(obj, X, Y, 128, 2);
            
            % TODO: more efficient mode collection
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
                    cx = cX(:, start:finish);
                    cy = cY(:, start:finish);
                    
                    if rank > winLen 
                        rank = winLen;
                    end
                    
                    [Ur, Sr, Vr, Atilde] = fitKoopman(obj, cx, cy, rank, tSize);
                    [lambda, omega, Phi, b] = fitModes(obj, Ur, Sr, Vr, Atilde, cx, cy, 50, 'exact-mod', x, y);
                    lambdaCum = [lambdaCum; lambda];
                    omegaCum = [omegaCum; omega];
                    PhiCum = [PhiCum, Phi];
                    bCum = [bCum; b];
                end
            end
            
            [lambdaCum, omegaCum, PhiCum, bCum, ~] = cleanModes(obj, lambdaCum, omegaCum, PhiCum, bCum);
            [norms, mags, freqs] = interpModes(obj, lambdaCum, PhiCum, bCum);
            
            [norms, mags, freqs, lambdaCum, omegaCum, PhiCum, bCum] = boxSelection(obj, norms, mags, freqs, ...
                                                                                    lambdaCum, omegaCum, PhiCum, bCum, ....
                                                                                    [0, 1], [0.9, 1], [0.8,15]);
            
            magFreq(obj, norms, mags, freqs);
            domStruct(obj, freqs, lambdaCum, PhiCum, bCum, tSize);
            nyqPlot(obj, norms, freqs, mags, lambdaCum);
            
            powSpect(obj, X, freqs, mags, 50);
        end
        
        % TODO: mode selection for EIT and improve mode selection for b
        function [norms, mags, freqs, lambda, omega, Phi, b] = boxSelection(obj, norms, mags, freqs, lambda, omega, Phi, b, limsNorms, limsMags, limsFreqs)
            nnorms = norms / max(norms);
            
            isNorms = nnorms >= limsNorms(1) & nnorms <= limsNorms(2);
            isMags = mags >= limsMags(1) & mags <= limsMags(2);
            isFreqs = freqs >= limsFreqs(1) & freqs <= limsFreqs(2);
            
            isSelect = isNorms' & isMags & isFreqs;
            
            norms = norms(isSelect);
            mags = mags(isSelect);
            freqs = freqs(isSelect);
            
            lambda = lambda(isSelect);
            omega = omega(isSelect);
            Phi = Phi(:, isSelect);
            b = b(isSelect);
        end
    end
end

