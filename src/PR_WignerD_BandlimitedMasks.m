% Copyright (c) 2018- Michigan State University, University of Michigan - 
% Dearborn and the CHARMS Research Group
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


%% (Ptychographic) Phase Retrieval using Wigner Deconvolution (Bandlimited Masks)
%
% Script to implement phase retrieval using Wigner Deconvolution for
% bandlimited masks
%
% For a description of the method, see
%
% 'Inverting Spectrogram Measurements via Aliased Wigner Distribution
% Deconvolution and Angular Synchronization' 
% Michael Perlmutter, Sami Merhi, Aditya Viswanthan and Mark Iwen
% https://arxiv.org/abs/1907.10773
%

clear; close all; clc

% For repeatability, set random seed to ...
rng(1234);


%% Parameters
% Signal parameters
d =  2^06;                      % signal dimension

% Mask parameters
% maskType = 'random';            % random masks
maskType = 'exp';               % (Non-symmetric) exponential mask from 
                                % earlier paper

rho = ceil(1.25*log2(d));     % support of m^hat
kappa = ceil(rho/4);          % Note: kappa should obey L=rho-1+kappa
L = rho+kappa-1;              % subsampling in space - no. of shifts

% Noise parameters
addnoise = false;               % Add noise?
snr = 50;                       % SNR of noise to be added


% Fix d to satisfy divisibility requirements
d = L*ceil(d/L);


% Print out problem parameters
fprintf( '\n\n --------------------------------------------------- \n' );
fprintf(     '    Phase Retrieval using Wigner Deconvolution \n' );
fprintf(     ' --------------------------------------------------- \n' );

fprintf( ' Problem size, d = %d \n', d );
fprintf( ' Bandwidth of mask (card. of support of m^hat), rho = %d \n', ...
                                     rho );
                                 
fprintf( ' No. of shifts (in space), L = %d \n', L );
fprintf( ' Total no. of measurements, D = %d = %dd\n\n', d*L, L );

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



%% Choice of signal x
x = randn(d,1)+1i*randn(d,1); % random complex signal x


%% Choice of mask m

% define m^hat
m_hat = zeros(d,1);

% Generate mask depending on type
switch lower(maskType)
    case 'random'
        % Random mask
        m_hat(1:rho) = randn(rho,1) .* exp(1i*2*pi*rand(rho,1));
    case 'exp'        
        % Exponential mask (deterministic, real)
        % parameter
        a = max(4, (rho-1)/2);
        m_hat(1:rho) = exp(-(0:rho-1).'/a)/((2*rho-1)^.25);
end
    
% here is the mask in physical space (for reference)
m = ifft(m_hat);


% Pre-computation (of terms involving masks)
mask_precomp_ld = zeros(d,kappa);
mask_precomp_ud = zeros(d,kappa-1);

for alpha = 0:kappa-1
    mask_precomp_ld(:,alpha+1) = fft( m_hat.*circshift(conj(m_hat),-alpha), d );
end

lb_idx = L-(kappa-1);
for alpha = L-(kappa-1):L-1
    mask_precomp_ud(:,alpha-lb_idx+1) = fft( m_hat.*circshift(conj(m_hat),L-alpha), d );
end


%% Measurements
Y = zeros(d,L);             % matrix of (spectrogram) measurements

% Each column   of Y corresponds of a physical shift \ell
% Each row      of Y corresponds to a Fourier  mode  \omega
% TODO: <change this to only compute subsampled measurements>
for l=0:d/L:d-1
    Y(:,l*(L/d)+1) = abs( fft(x.*circshift(m,l), d) ).^2;
end

% TODO:
% Y_alt = spectrogram(x, m, rho-1, d);

% noise
if(addnoise)
    signal_power = norm( Y(:) )^2 / (d*L);
    noise_power = signal_power / ( 10^(snr/10) );    
    % Add (real) Gaussian noise of desired variance
    noise = sqrt(noise_power)*randn( size(Y) );
    Y = Y + noise;    
else
    noise_power = 0;
end


%% Solve for diagonal of x^x^*
tic;
% First, compute the left-hand side of double-aliasing formulation
% TODO: add Sami's paper reference
LHS = fft( fft(Y,L,2),d ).';

% Next, solve for diagonals
X_diags = zeros(d,2*kappa-1);
for alpha = 0:kappa-1
    X_diags(:,alpha+1) =  LHS(alpha+1,:).' ./ ...
                                            mask_precomp_ld(:,alpha+1);  
end

lb_idx = L-(kappa-1);
for alpha = L-(kappa-1):L-1    
    X_diags(:,kappa+alpha-lb_idx+1) = LHS(alpha+1,:).' ...
                                    ./ mask_precomp_ud(:,alpha-lb_idx+1); 
end


% Todo: fix the unnecessary conj(.)
X_diags = conj(ifft(X_diags, d));


% Put recovered entries into diagonal form
diagLocs = [0:kappa-1 -(kappa-1):-1 ...
            -(d-1):-(d-1)+(kappa-2) (d-1)-(kappa-2):(d-1)];
X = spdiags([X_diags X_diags(:,2:end)], diagLocs, d, d);


% Hermitian symmetrize
X = X/2 + X'/2;


% Note down magnitudes
mags = (d/sqrt(L))*sqrt( diag(X) );


% entrywise normalize to get a matrix of relative phases
nz_idx = find(X);
X(nz_idx) = X(nz_idx)./abs(X(nz_idx));


%% Angular synchronization

% this is Matlab's eig solver
[xrec, ~, ~] = eigs(X, 1, 'LM');    % compute leading eigenvector
xrec = xrec./abs(xrec);             % retain only phase

% we have recovered the DFT of x; now get back x
xrec = ifft( full(mags.*xrec), d );


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

    