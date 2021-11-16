classdef DMD < handle
    %DMD Summary of this class goes here
    %   Detailed explanation goes here
    
    % TODO
%     Implement MrDMD
%     What is b?
%     Why scale lambda with log/sub/dt/2/pi?
%     Try running DMD at each sample
%     How long is a prediction available for?
%     Define lengths of windows
    
    properties
        dt
        Fs
        t
        X2
        X1
        Rank
        Atilde
        Ur
        Sr
        Vr
        W
        Lambda
        Phi
        Omega
        B
        Xdmd
        Dynamics
    end
    
    methods
        function obj = DMD()
            %DMD Construct an instance of this class
        end
        
        function loadData(obj,X1,X2,fs)
            obj.X1 = X1;
            obj.X2 = X2;
            obj.Fs = fs;
            obj.dt = 1/fs;
            obj.t = (0:size(X1,2) - 1) * obj.dt;
        end
        
        function truncate(obj,rank) % TODO: Optimal SVD truncation for EIT
            if isempty(rank)
                [U,S,V] = svd(obj.X1,'econ');
                beta = size(obj.X1,2)/size(obj.X1,1);
                sigma = diag(S);
                tau = optimal_SVHT_coef(beta,0) * median(sigma);
                
                rank = find(sigma > tau,1,'last');
                obj.Ur = U(:,1:rank);
                obj.Sr = S(1:rank,1:rank);
                obj.Vr = V(:,1:rank);
                obj.Atilde = obj.Ur'*obj.X2*obj.Vr/obj.Sr;
            elseif rank > 0
                [U,S,V] = svd(obj.X1,'econ');
                obj.Ur = U(:,1:rank);
                obj.Sr = S(1:rank,1:rank);
                obj.Vr = V(:,1:rank);
                obj.Atilde = obj.Ur'*obj.X2*obj.Vr/obj.Sr;
            else
                obj.Atilde = obj.X2\obj.X1;
            end
            obj.Rank = rank;
        end
        
        function [] = fit(obj) % TODO: Total-least-quares or forward/backward
            [W,lambda] = eig(obj.Atilde);
            
            obj.Lambda = diag(lambda);
            obj.Omega = log(obj.Lambda)/obj.dt;
            
            obj.Phi = obj.X2*obj.Vr/obj.Sr*W;
            
            obj.B = obj.Phi\obj.X1(:,1);
            obj.W = W;
            % or
            % alpha1 = S(1:r,1:r)*V(1,1:r)’;
            % bPOD = ( Atilde*W)\alpha1;
        end
        
        function getDynamics(obj)
            dynamics = zeros(obj.Rank,length(obj.t));
            
            for i = 1:length(obj.t)
                dynamics(:,i) = (obj.B .* exp(obj.Omega * obj.t(i)));
            end
            
            obj.Dynamics = dynamics;
        end
        
        function dmdSpectrum(obj)
            Ahat = (obj.Sr ^ (-1/2)) * obj.Atilde * (obj.Sr ^ (1/2));
            [What, ~] = eig(Ahat);
            Wr = obj.Sr ^ (1/2) * What;
            phi = obj.X2 * obj.Vr / obj.Sr * Wr;
            
            f = abs(imag(obj.Omega));
            P = diag(phi' * phi);
            
            figure;
            stem(f, P, 'k');
            xlim([0 obj.Fs/2]);
        end
        
%         function dmdSpectrumVandermonde(obj)
%             Vand = zeros(obj.Rank, size(obj.X1, 2)); % Vandermonde matrix
%             for k = 1: size(obj.X1, 2)
%                 Vand(:, k) = obj.Lambda.^(k-1);
%             end
%             
%             G = obj.Sr * obj.Vr';
%             P = ( obj.W'*obj.W).*conj(Vand*Vand');
%             q = conj(diag(Vand*G'*obj.W));
%             Pl = chol(P,'lower');
%             b = ( Pl')\(Pl\q);
%             
%             f = abs(imag(obj.Omega));
%             figure;
%             stem(f, b, 'k');
%             xlim([0 obj.Fs/2]);
%         end
        
        function powerSpectrum(obj)
            timesteps = size(obj.X1, 2); 
            srate = obj.Fs;
            nelectrodes = size(obj.X1,1);
            NFFT = 2^ nextpow2(timesteps);
            f = srate/2*linspace(0, 1, NFFT/2+1);
            
            figure
            hold on;
            for c = 1: nelectrodes
                fftp(c,:) = fft(obj.X1(c,:), NFFT);
                plot(f, 2* abs(fftp(c,1:NFFT/2+1)), ...
                    'Color', 0.6*[1 1 1]);
            end
%             plot(f, 2* abs(mean(fftp(c,1:NFFT/2+1), 1)), ...
%                 'k', 'LineWidth', 2); 
            xlim([0 obj.Fs/2]);
%             ylim([0 400]);
        end
        
        function dmdPoles(obj)
            figure;
            subplot(1,2,1)
            plot(obj.Lambda, 'k.');
            rectangle('Position', [-1 -1 2 2], 'Curvature', 1, ...
                'EdgeColor', 'k', 'LineStyle', '--');
            
            subplot(1,2,2)
            plot(obj.Omega, 'k.');
            line([0 0], 200*[-1 1], 'Color', 'k', 'LineStyle', '--');
        end
        
        function mrFit(obj, cutoff, nLevels)
            for i = 1:nLevels
                fit(obj);
%                 removeIdx = 
                filterModes(obj,idx);
            end
        end
        
        % TODO: mode selection for EIT
        function selectModes(obj)
            
        end
        
        function Xdmd = filterModes(obj,idx)
            phiFilt = obj.Phi;
            phiFilt(:,idx) = 0;
            
            Xdmd = phiFilt * obj.Dynamics;
        end
    end
end

