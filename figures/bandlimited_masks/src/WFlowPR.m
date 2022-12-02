function [xrec, etime] = WFlowPR(Y, pars, mask)

%% Initialization
% Read in problem parameters
d = pars.d;
L = pars.L;
m = mask.m;

snr = pars.snr;
               

%% Phase Retrieval using Wirtinger Flow

% Masks for WF
Masks = zeros(d,L);
for ix = 1:L
    Masks(:,ix) = circshift(m, (ix-1)*d/L );
end

% Make linear operators; A is forward map and At its scaled adjoint (At(Y)*numel(Y) is the adjoint) 
A_CDP = @(I)  fft(Masks .* repmat(I,[1 L]));  % Input is n x 1 signal, output is n x L array
At_CDP = @(Y) mean(conj(Masks) .* ifft(Y), 2);            % Input is n x L array, output is n x 1 signal          

    
%% Wirtinger Flow -- CDP 

tic;

% Initialization
npower_iter = 50;                           % Number of power iterations 
z0 = randn(d,1); z0 = z0/norm(z0,'fro');    % Initial guess 
for tt = 1:npower_iter 
    z0 = At_CDP(Y.*A_CDP(z0)); z0 = z0/norm(z0,'fro');
end

normest = sqrt(sum(Y(:))/(d*L));    % Estimate norm to scale eigenvector  
z = normest * z0;                   % Apply scaling 


%% Loop

% Don't know how to choose these - seems to vary with d and other
% parameters
switch snr
    case 10
        T = 150;
    case 20
        T = 500;
    case 30
        T = 1200;
    case 40
        T = 2000;
    case 50
        T = 3000;
    case 60
        T = 4500;
end


% T = 2500;                           % Max number of iterations
tau0 = 300;                         % Time constant for step size
mu = @(t) min(1-exp(-t/tau0), 0.2);

for t = 1:T
    Bz = A_CDP(z);
    C  = (abs(Bz).^2-Y) .* Bz;
    grad = At_CDP(C)*d;                    % Wirtinger gradient
    z = z - mu(t)/normest^2 * grad;  % Gradient update 
end

etime = toc;
xrec = z;

return