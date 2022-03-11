addpath ./dmdCode/CH09_SPARSITY/utils
% PARAMETERS
n = 128;
K = 5; % degree of sparsity
T = 2; % duration of integration
dt = 0.01; % time step
r = 10; % number of PCA modes to keep
noisemag = 0.0; % magnitude of additive white noise

xtilde = zeros(n,n);
% generate sparse FFT data
for i=1:K;
loopbreak = 0;
while(~loopbreak)
I(i) = 0+ceil(rand (1)*n/15) ;
J(i) = 0+ceil(rand (1)*n/15) ;
if(xtilde(I(i),J(i)) == 0)
loopbreak = 1;
end
end
IC(i) = randn();
xtilde(I(i),J(i)) = IC(i);
F(i) = sqrt (4*rand());
damping(i) = -rand ()*.1;
end

xtilde = zeros(n,n);
XDAT = [];
XDATtilde = [];
XDATtildeNoise = [];
for t = 0:dt:T
xtilde = 0*xtilde;
for k=1:K
xtilde(I(k),J(k)) = exp(damping(k)*t)*cos(2*pi*F(k)*t)*...
IC(k) + exp(damping(k)*t)*sqrt (-1)*sin(2*pi*F(k)*t)*...
IC(k);
end
XDATtilde = [XDATtilde reshape(xtilde ,n^2,1)];
xRMS = sqrt ((1/(n*n))*sum(xtilde(:).^2));
xtilde = xtilde + noisemag*xRMS *randn(size(xtilde)) +...
noisemag*xRMS*sqrt (-1)*randn(size(xtilde));
XDATtildeNoise = [XDATtildeNoise reshape(xtilde ,n^2,1)];
x = real (ifft2(xtilde));
XDAT = [XDAT reshape(x,n^2,1)];
end

M = 15; % number of measurements
% projType = 1; % uniform random projection
projType = 2; % Gaussian random projection
% projType = 3; % Single pixel measurement
%% Create random measurement matrix
C = zeros(M,n*n);
Theta = zeros(M,n*n);
xmeas= zeros(n,n);
for i=1:M
xmeas = 0*xmeas;
if(projType==1)
xmeas = rand(n,n);
elseif(projType==2)
xmeas = randn(n,n);
elseif(projType==3)
xmeas(ceil (n*rand),ceil(n*rand)) = 1;
end
C(i,:) = reshape(xmeas ,n*n,1);
Theta(i,:) = reshape((ifft2(xmeas)),1,n*n);
end

%% Project data
YDAT = C*XDAT;
%% plot C
figure , colormap bone
Cimg = sum(C,1);
Ximg = XDAT (:,1).*Cimg';
imagesc(reshape(Ximg,n,n));
figure , colormap bone;
surf(X,Y,Z,reshape(Cimg ,n,n));
view (10,32)
axis equal , axis off
drawnow
