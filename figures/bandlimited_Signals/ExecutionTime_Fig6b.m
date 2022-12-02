% Copyright (c) 2018- Michigan State University and the CHARMS Research
% Group
% 
% This file is part of the BlockPR software package
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy of
% this software and associated documentation files (the "Software"), to deal in
% the Software without restriction, including without limitation the rights to
% use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
% of the Software, and to permit persons to whom the Software is furnished to do
% so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

%% (Ptychographic) Phase Retrieval using Wigner Deconvolution
%
% Script to implement phase retrieval using using Wigner Deconvolution. 
%   
% For a description of the method, see arXiv preprint
% https://arxiv.org/abs/1907.10773
%

% clear; close all; clc
clear;
addpath src/

% For repeatability, set random seed to ...
rng(12345);


%% Parameters

% Use the following values of d, gamma
% Note: Cannot choose d arbitrarily - has to satisfy divisibility
% conditions; moreover, certain values of (d,gamma) allow for better
% conditioned Vandermonde systems
dvals = [38; 114; 266; 418; 570; 722; 874; 1026];
gammavals = [10; 29; 67; 105; 143; 181; 219; 257];

% Store execution times here
etime = zeros( length(dvals), 3 );

% No. of trials
ntrials = 100;


% For each value of d...
for idval = 1:length(dvals)    
    
d =  dvals(idval);                 % signal dimension
gamma = gammavals(idval);          % support of xhat
delta = 10;                         % support of m

addnoise =  false;              % Add noise?
snr = 50;                       % SNR of noise to be added

kappa = round(delta/1);         % Note: kappa should obey L=delta-1+kappa
xi = round(gamma/1);
K = delta+kappa-1;              % subsampling in frequencies - no. of modes
L = gamma+xi-1;                 % subsampling in space - no. of shifts

fprintf( '\n Now running simulations for d = %d ...\n', d );


%% Choice of mask m
% Choose "best" random mask (based on mask constant)

% initialization
idx = 1; mu = 0;

while( idx<=100 || mu<=1e-5 )
    % define m^hat
    m = zeros(d,1);
    m(1:delta) = (1+0.5*rand(delta,1)) .* exp(1i*2*pi*rand(delta,1));    

    % here is the mask in physical space (for reference)
    mhat = fft(m);

    % Pre-computation (of terms involving masks)
    mask_precomp_ld = zeros(2*delta-1,gamma);
    mask_precomp_ud = zeros(2*delta-1,gamma-1);

    for alpha = 0:gamma-1
        tmp = fft( mhat.*circshift(conj(mhat),-alpha), d );
        tmp = [tmp(1:delta); tmp(d-delta+2:end)];
        mask_precomp_ld(:,alpha+1) = tmp(1:end);
    end

    for alpha = gamma:2*gamma-2
        tmp = fft( mhat.*circshift(conj(mhat),L-alpha), d );
        tmp = [tmp(1:delta); tmp(d-delta+2:end)];
        mask_precomp_ud(:,alpha-gamma+1) = tmp(1:end);
    end
    
    % mu value
    mu_current = min(abs([mask_precomp_ld(:); mask_precomp_ud(:)]));
    
    % update best choice of masks
    if(mu<mu_current)
        mu = mu_current;
        mhat_best = mhat;
        mbest = m;
        mask_precomp_ld_best = mask_precomp_ld;
        mask_precomp_ud_best = mask_precomp_ud;
    end
        
    % update loop counter
    idx = idx+1;
end

% % Display mask-constant mu
% mu

% Store these in a structure
mask = struct(  'ld', mask_precomp_ld_best, ...
                'ud', mask_precomp_ud_best, ...
                 'm', mbest );

% use the "best" mask
m = mbest;


% Fourier transform of mask
mhat = fft(m);


% Vandermonde matrix (see Step 4 in Alg. 2)
W = zeros(2*delta-1,gamma);
for j = 0:2*delta-2
    for k = 0:gamma-1
        W(j+1,k+1) = exp(-2i*pi*(j-delta+1)*k/d);
    end
end
% cond(W)


% for the HIO+ER implementation (note: could be made more efficient) 
% Sub-sampled Fourier matrix
Fd = dftmtx(d);
Fk = Fd(1:d/K:end,:);
clear Fd;

% Measurement matrix
M = zeros(L*K, d);
ctr = 1;
for ix = 1:d/L:d
    M((ctr-1)*K+1:ctr*K, :) = Fk*diag(circshift(m,ix-1));
    ctr = ctr+1;
end
Mproj = M*pinv(M);


% Reg. parameter for the iterated Tikhonov procedure (changes with 
% the noise level; use L-curve to find this)
% ... for Alg. 2 only
regpar_alg2 = 1e2;
% ... for Alg. 2 + Alg. 3
regpar = 1e2;


% Initialize regularization and iteration parameters
alpha0 = regpar;
q = 8/10;

    
% How do we invert the Vandermonde system?
% For Alg. 2
pinv_W = (W'*W + regpar_alg2*eye(gamma))\W';

% For Alg. 2 (w/ steps (4)-(6) replaced by Alg. 3)
pinvW_Tikhonov = zeros(gamma, 2*delta-1, 21);

for ctr = 0:20
    pinvW_Tikhonov(:,:,ctr+1) = (W'*W + (alpha0*q^ctr)*eye(gamma))\W';
end

% For inverting the Fourier system (note: efficiency can be improved using
% FFTs)
F_L = dftmtx(L);
F_K = dftmtx(K);

% For circular indexing
wrap_d = @(x) (1 + mod(x, d));
wrap_K = @(x) (1 + mod(x, K));
wrap_L = @(x) (1 + mod(x, L));


%% For each trial...
for itrial = 1:ntrials


%% Choice of signal x
% (bandlimited) random complex signal x 
xhat = zeros(d,1);
xhat(1:gamma) = randn(gamma,1)+1i*randn(gamma,1);
x = ifft(xhat);


%% Measurements
Y = zeros(K,L);             % matrix of (spectrogram) measurements

% Each column   of Y corresponds of a physical shift \ell
% Each row      of Y corresponds to a Fourier  mode  \omega
% TODO: <change this to only compute subsampled measurements>

for l=0:d/L:d-1
    tmp = abs( fft(x.*circshift(m,l), d) ).^2;
    Y(:,l*(L/d)+1) = tmp(1:d/K:end);
end

% TODO:
% do this using spectrogram command

% noise
if(addnoise)
    signal_power = norm( Y(:) )^2 / (L*K);
    noise_power = signal_power / ( 10^(snr/10) );    
    % Add (real) Gaussian noise of desired variance
    noise = sqrt(noise_power)*randn( size(Y) );
    Y = Y + noise;    
else
    noise_power = 0;
end


%% Solve for diagonal of xx*
tic;

% First, compute the left-hand  side of double-aliasing formulation

%%%%%% To do: Efficiency can be improved using FFTs
% LHS = fft( fft(Y,L,2),K ).';
LHS = F_L*Y.'*F_K.';


% This corresponds to Steps 2 and 3 in Alg. 2
V = zeros(2*delta-1, 2*gamma-1);


for omega = 0:delta-1
    for alpha = 0:gamma-1
        
        tmp = fft( mhat .* conj( circshift(mhat, -alpha) ) );
        
        RHS = (d^3/(K*L)) * LHS( wrap_L(alpha), wrap_K(omega) ) / ...
                                        tmp( wrap_d(omega) );  
                                    
        V( delta+omega, gamma-alpha ) = RHS;
    end
    
    
    for alpha = gamma:2*gamma-2
        
        tmp = fft( mhat .* conj( circshift(mhat, L-alpha) ) );
        
        RHS = (d^3/(K*L)) * LHS( wrap_L(alpha), wrap_K(omega) ) / ...
                                        tmp( wrap_d(omega) );  
                                    
        V( delta+omega, 3*gamma-(alpha+1) ) = RHS;
    end    
    
end
        


for omega = delta:2*delta-2
    for alpha = 0:gamma-1
        
        tmp = fft( mhat .* conj( circshift(mhat, -alpha) ) );
        
        RHS = (d^3/(K*L)) * LHS( wrap_L(alpha), wrap_K(omega) ) / ...
                                        tmp( wrap_d(omega-K) );  
                                    
        V( omega-delta+1, gamma-alpha ) = RHS;
    end
    
    
    for alpha = gamma:2*gamma-2
        
        tmp = fft( mhat .* conj( circshift(mhat, L-alpha) ) );
        
        RHS = (d^3/(K*L)) * LHS( wrap_L(alpha), wrap_K(omega) ) / ...
                                        tmp( wrap_d(omega-K) );  
                                    
        V( omega-delta+1, 3*gamma-(alpha+1) ) = RHS;
    end    
    
end

% Since we will use different methods of inverting the Vandermonde system, 
% we keep track of this inversion times separately
part1_etime = toc;
tic; 


%% Iterated Tikhanov regularization (see Alg. 3)
% Some temporary variables
A_skew = zeros(gamma,2*gamma-1);
G = zeros(gamma);

resdl = ones(size(V));
ctr = 0;
alpha0 = regpar;

% We do 20 steps of the iteratoin
while( ctr<=20 )    

    % First, consider a rank-1 projection
    [Ug,Sg,Vg] = svd(G);
    Sg(2:end) = 0;
    A = Ug*Sg*Vg';

    % Convert b/w column and diagonal forms
    A_skew = zeros(gamma,2*gamma-1);
    for ix = 1:gamma
        A_skew(ix,gamma-ix+1:2*gamma-ix) = A(ix,:);
    end

    % Residual after inverting Vandermonde system
    resdl = V - W*A_skew;
    % Apply regularization
    A_skew = A_skew + pinvW_Tikhonov(:,:,ctr+1)*resdl;
    ctr = ctr+1;

    % Convert b/w column and diagonal forms    
    A = zeros(gamma);
    for ix = 1:gamma
        A(ix,:) = A_skew(ix, gamma-ix+1:2*gamma-ix);
    end

    % Hermitian symmetrize
    G = A/2 + A'/2;

end

%% Angular synchronization
% this is Matlab's eig solver
% Note: the matrix G is typically small; if not, eigs or other alternatives
% could be used
[evec, eval] = eig(G);
[lambdamax, loc] = max(diag(abs(eval)));

xhat_e = zeros(d,1);
xhat_e(1:gamma) = sqrt(lambdamax)*evec(:,loc);
x_e = ifft(xhat_e);


% Correct for global phase factor
phaseOffset = angle( (x_e'*x) / (x'*x) );
xrec = x_e * exp(1i*phaseOffset);

% Store execution time
etime(idval, 1) = etime(idval, 1) + (part1_etime+toc)/ntrials;


%% This is for the naive method with no fancy regularization
    tic;
    % Invert Vandermonde system 
    A_skew = pinv_W*V;

    % Convert b/w column and diagonal forms    
    A = zeros(gamma);
    for ix = 1:gamma
        A(ix,:) = A_skew(ix, gamma-ix+1:2*gamma-ix);
    end

    % Hermitian symmetrize
    G = A/2 + A'/2;


    %% Angular synchronization
    % this is Matlab's eig solver
    [evec, eval] = eig(G);
    [lambdamax, loc] = max(diag(abs(eval)));

    xhat_e = zeros(d,1);
    xhat_e(1:gamma) = sqrt(lambdamax)*evec(:,loc);
    x_e = ifft(xhat_e);

    % Correct for global phase factor
    phaseOffset = angle( (x_e'*x) / (x'*x) );
    xrec = x_e * exp(1i*phaseOffset);

    % Store execution time
    etime(idval, 2) = etime(idval, 2) + (part1_etime+toc)/ntrials;



%% HIO+ER Reconstruction
tic; 
xrec = postprocess_HIOER( zeros(d,1), Y, m, d, K, L, [600, 25, 5], M, Mproj );

% Correct for global phase factor
phaseOffset = angle( (xrec'*x) / (x'*x) );
xrec = xrec * exp(1i*phaseOffset);

% Store execution time
etime(idval, 3) = etime(idval, 3) + toc/ntrials;


end

fprintf( '\n  (Avg.) execution time for WignerD (Alg. 2 + Alg. 3) is %3.2es', ...
                                    etime(idval, 1) );                 
fprintf( '\n  (Avg.) execution time for WignerD (Alg. 2) is %3.2es', ...
                                    etime(idval, 2) );                                
fprintf( '\n  (Avg.) execution time for HIO+ER is %3.2es\n', ...
                                    etime(idval, 3) );

end


% plot
% Note: you may need to change the axis limits depending on the speed of 
% your hardware/software implementation
figure; 
loglog(dvals, etime(:,1), '-+', 'linewidth', 2); hold on
loglog(dvals, etime(:,2), '--o', 'linewidth', 2);
loglog(dvals, etime(:,3), '-.s', 'linewidth', 2);
axis([1e1 2e3 1e-3 1e2]); grid
xlabel( 'Problem Size, $d$', 'interpreter', 'latex', 'fontsize', 14 );
ylabel( 'Execution Time (secs.)', 'interpreter', 'latex', 'fontsize', 14 );
title( 'Execution Time ($\gamma\approx d/4,\delta=10$)', ...
                                'interpreter', 'latex', 'fontsize', 14 );
legend( {'Alg. 2 (w/ Steps (4)-(6) replaced by Alg. 3)', 'Alg. 2', ...
                    'HIO+ER'}, ...
         'location', 'northwest', 'interpreter', 'latex', 'fontsize', 14 );


dt = linspace(1.25e2,5.5e2,50);                             
loglog( dt, 3e-5*(dt.^2), 'k--', 'linewidth', 2 ,'HandleVisibility','off'); 
loglog( dt, 5e-6*dt.*log2(dt), 'k--', 'linewidth', 2 ...
                                    ,'HandleVisibility','off');
drawnow;
                                
                
% Save data
% save( 'etime.mat', 'dvals', 'etime' );

fprintf( '\n\n ' );
                                 