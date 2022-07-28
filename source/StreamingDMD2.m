%StreamingDMD2 
%Liew

classdef StreamingDMD2 < handle
    %calculate DMD in streaming mode
    properties
        rmin;
        rmax;
        thres;
        halflife;
        rho;
        Ux;
        Uy;
        Pinvx;
        Pinvy;
        Q;
    end   
        
    methods

        
        
        
        function obj = StreamingDMD2(X, Y,rmin, rmax,thres, halflife)
            
            obj.rmin = rmin;
            obj.rmax = rmax;
            obj.thres = thres;
         
            obj.halflife = halflife;
            if halflife == 0
                obj.rho = 1;
               
            else obj.rho = 2*(-1/halflife);
            end
            
            %Eq(2) truncated SVD
%             [obj.Ux,~,~] = truncSVD(X, obj.rmin);
%             [obj.Uy,~,~] = truncSVD(Y, obj.rmin);

            [Ux,Sx,Vx] = svd(X,'econ');
            obj.Ux = Ux(:,1:obj.rmin);
            Srx = Sx(1:obj.rmin,1:obj.rmin);
            Vrx = Vx(:,1:obj.rmin);
            
            [Uy,Sy,Vy] = svd(Y,'econ');
            obj.Uy = Uy(:,1:obj.rmin);
            Sry = Sy(1:obj.rmin,1:obj.rmin);
            Vry = Vy(:,1:obj.rmin);
            
            %Eq(10) mapping of input vector to reduced order space
            X_tild = obj.Ux.' * X;
            
            %Eq(10) mapping of out vector to reduced order space
            Y_tild = obj.Uy.' * Y;
            
            %Eq(9) Decomposition of transition matrix into product of Q and
            %Pinv
            obj.Q = Y_tild * X_tild.';
            obj.Pinvx = X_tild * Y_tild.';
            obj.Pinvy = Y_tild * Y_tild.';
        end 
        
        function update(obj, x, y)
            % RESHAPE INTO COLUMNS
%              x, y = x.reshape([-1, 1]), y.reshape([-1, 1])
%              status = 0
% 
            % GET 2-NORMS
            normx = norm(x);
            normy = norm(y);
%           
            % PROJECT NEW VECTORS ONTO SINGULAR VECTORS OF LARGE MATRICES
            xtilde = obj.Ux.' * x; % Ux' * x
            ytilde = obj.Uy.' * y; % Uy' * y
            
            % CALCULATE ERROR OF THE PROJECTION
            %Numeraor of Eq.(14) projection error
            ex = x - obj.Ux * xtilde;
            ey = y - obj.Uy * ytilde;
            
%             x_status = Status.NONE
%         y_status = Status.NONE

            %%%STEP 1 BASIS EXPANSION%%%
            % INCREASE RANK OF Ux AND Uy WHEN PROJECTION ERROR RISES
            %Table 1: Rank augmentation of Ux
            
            if norm(ex)/normx > obj.thres
                % GET NEW VECTOR FOR Ux MATRIX AND APPEND TO SECOND
                % DIMENSION
                u_new = ex / norm(ex);
                u_new_r = reshape(u_new,[],1);
                obj.Ux = cat(2,obj.Ux,u_new_r);
                
                obj.Pinvx = cat(2,obj.Pinvx,zeros(size(obj.Pinvx,1),1));
                obj.Pinvx = cat(1,obj.Pinvx,zeros(1,size(obj.Pinvx,2)));
                
                %obj.Q = cat(1,zeros(1,size(obj.Q,2)));
                obj.Q = cat(2,zeros(size(obj.Q,1),1));
               % x_status = Status.AUGMENTATION; %?????
            end
            
            %Table 1: Rank augmentation of Uy
            
            if norm(ey)/normy > obj.thres
                u_new = ey / norm(ey);
                u_new_r = reshape(u_new,[],1);
                obj.Uy = cat(2,obj.Uy,u_new_r);
                
                obj.Pinvy = cat(2,obj.Pinvy,zeros(size(obj.Pinvy,1),1));
                obj.Pinvy = cat(1,obj.Pinvy,zeros(1,size(obj.Pinvy,2)));
                
                %obj.Q = cat(2,zeros(size(obj.Q,1),1));
                obj.Q = cat(1,obj.Q, zeros(1,size(obj.Q,2)));
                
               % y_status = Status.AUGMENTATION; %?????
            end
            
            %%%StEP 2 - BASIS POD COMPRESSION
            %Table 1: Rank reduction of Ux
            
            if rank(obj.Ux) > obj.rmax
                [eigvec, D] = eig(obj.Pinvx);
                eigval = diag(D);
                [~,indx] = sort(-eigval);
                eigval = eigval(indx);
                qx = eigvec(:,indx(1:obj.rmin));
                
                obj.Ux = obj.Ux * qx;
                obj.Q = obj.Q * qx;
                obj.Pinvx = diag(eigval(1:obj.rmin));
                %x_status = Status.REDUCTION;
            end
            
            %Table 1: Rank reduction of Uy
            if rank(obj.Uy) > obj.rmax
                [eigvec, D] = eig(obj.Pinvy);
                eigval = diag(D);
                [sorted,indx] = sort(-eigval);
                eigval = -sorted;
                qy = eigvec(:,indx(1:obj.rmin));
                
                obj.Uy = obj.Uy * qy;
                obj.Q = qy.' * obj.Q;
                obj.Pinvy = diag(eigval(1:obj.rmin));
                %x_status = Status.REDUCTION;
            end
           
            
            %%%STEP 3 REGRESSION UPDATE
            xtilde = obj.Ux.' * x;
            ytilde = obj.Uy.' * y;
            
            %Eq 10, 11 and 12 - Rank 1 update of DMD matrices
            obj.Q = obj.rho * obj.Q + ytilde * xtilde.';
            obj.Pinvx = obj.rho * obj.Pinvx + xtilde * xtilde.';
            obj.Pinvy = obj.rho * obj.Pinvy + ytilde * ytilde.';
            
            %return x_status, y_status
        end
        
        %%@property
        function Ushape = rank(obj)
            Ushape = size(obj.U.shape);
        end
        
        function temp = A(obj)
            temp = obj.Q * pinv(obj.Pinvx);
        end
        
%         function [modes,eigvals] = modes(obj)
%             [eigvals, eigvecK] = eig(obj.Ux.' * obj.Uy * obj.A);
%             modes = obj.Ux * eigvecK;
%         end
        
        function[modes,eigvals] = modes(obj)
            [V,D] = eig(obj.Ux.' * obj.Uy * obj.A);
            modes = obj.Ux * V;
            eigvals = diag(D);
        end       
   
    end
end


  
            
                            
    
    