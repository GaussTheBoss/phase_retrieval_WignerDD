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

clear; close all; clc

% For repeatability, set random seed to ...
rng(12345);


%% Parameters

% Signal parameters
d =  114;                       % signal size
gamma = 29;                     % support of xhat (Signal is bandlimited)

% Mask parameters
% maskType = 'random';            % random masks
maskType = 'exp';               % (Non-symmetric) exponential mask from 
                                % earlier paper
delta = 10;                     % support of m 
kappa = round(delta/1);         % Note: kappa should obey L=delta-1+kappa
xi = round(gamma/1);
K = delta+kappa-1;              % subsampling in frequencies - no. of modes
L = gamma+xi-1;                 % subsampling in space - no. of shifts


% Noise parameters
addnoise = false;               % Add noise?
snr = 40;                       % SNR of noise to be added

% Fix d to satisfy divisibility requirements
d1 = K*round(d/K);
d2 = L*round(d/L);
d = d1;     % need to come up with a better procedure here


% Print out problem parameters
fprintf( '\n\n --------------------------------------------------- \n' );
fprintf(     '    Phase Retrieval using Wigner Deconvolution \n' );
fprintf(     ' --------------------------------------------------- \n' );

fprintf( ' Problem size, d = %d \n', d );
fprintf( ' Bandwidth parameter (for the bandlimited signal), gamma = %d \n', gamma );
fprintf( ' Spatial support of mask, delta = %d \n', delta );
fprintf( ' No. of Fourier modes (frequencies), K = %d \n', K );
fprintf( ' Total no. of measurements, D = %d = %dd\n\n', d*K, K );

switch lower(maskType)
    case 'random'
        fprintf( ' Using random masks \n\n' );
    case 'exp'
        fprintf( ' Using (non-symmetric) exponential masks from BlockPR method \n\n' );
end

if( addnoise )
    fprintf( ' Noisy simulation? - Yes \n' );
    fprintf( ' Added noise (SNR, dB) = %3.2f \n\n', snr );
else
    fprintf( ' Noisy simulation? - No \n\n' );
end


%% Choice of mask m
% initialize m
m = zeros(d,1);

% Generate mask depending on type
switch lower(maskType)
    case 'random'
        % Random mask
        m(1:delta) = randn(delta,1) .* exp(1i*2*pi*rand(delta,1));
    case 'exp'        
        % Exponential mask (deterministic, real)
        % parameter
        a = max(4, (delta-1)/2);
        m(1:delta) = exp(-(0:delta-1).'/a)/((2*delta-1)^.25);
end

% here is the Fourier transform of the mask
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
    

%% Other initialization (for the iterated Tikhonov regularization)
% Vandermonde matrix (see Step 4 in Alg. 2)
W = zeros(2*delta-1,gamma);
for j = 0:2*delta-2
    for k = 0:gamma-1
        W(j+1,k+1) = exp(-2i*pi*(j-delta+1)*k/d);
    end
end

% Some temporary variables
A_skew = zeros(gamma,2*gamma-1);
G = zeros(gamma);

% Initialize regularization and iteration parameters
% Reg. parameter for the iterated Tikhonov procedure (changes with 
% the noise level; use L-curve to find this)
alpha0 = 1e-10;     % Noiseless 
% alpha0 = 3e-0;         % Noisy (@40dB SNR)
q = 8/10;
    
% How do we invert the Vandermonde system? (for Alg. 3 - iterated Tikhonov)
pinvW_Tikhonov = zeros(gamma, 2*delta-1, 21);
for ctr = 0:20
    pinvW_Tikhonov(:,:,ctr+1) = (W'*W + (alpha0*q^ctr)*eye(gamma))\W';
end

% For circular indexing
wrap_d = @(x) (1 + mod(x, d));
wrap_K = @(x) (1 + mod(x, K));
wrap_L = @(x) (1 + mod(x, L));


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
LHS = fft( fft(Y,L,2),K ).';

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


%% Iterated Tikhanov regularization (see Alg. 3)
resdl = ones(size(V));
ctr = 0;

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
xrec = ifft(xhat_e);

% Correct for global phase factor
phaseOffset = angle( (xrec'*x) / (x'*x) );
xrec = xrec * exp(1i*phaseOffset);

etime = toc;
fprintf( ' Execution time is %3.3e secs.\n', etime );
                                    

%% Plot Reconstruction

% Physical space nodes
nodes = (0:d-1).';

% Plot the real part of signal
figure; subplot(1,2,1)
stem( nodes, real(x), 'ro', 'linewidth', 2 ); hold on
stem( nodes, real(xrec), 'b--+', 'linewidth', 2 );
xlabel n; ylabel x(n); 
title 'Phase Retrieval using Wigner Deconv. - Real Part'
legend('True', 'Recovered'); grid

% ... and the imaginary part
subplot(1,2,2)
stem( nodes, imag(x), 'ro', 'linewidth', 2 ); hold on
stem( nodes, imag(xrec), 'b--+', 'linewidth', 2 );
xlabel n; ylabel x(n); 
title 'Phase Retrieval using Wigner Deconv. - Imaginary Part'
legend('True', 'Recovered'); grid

% Reconstruction error
errordB = 10*log10( norm(x-xrec)^2 ...
                                / norm(x)^2 );
fprintf( '\n 2 Norm Error in reconstruction is %3.3e or %3.2f dB', ...
        norm(x-xrec)/norm(x), errordB );
fprintf( '\n Inf. Norm Error in reconstruction is %3.3e\n\n', ...
        norm(x-xrec,inf)/norm(x,inf) );

                                 