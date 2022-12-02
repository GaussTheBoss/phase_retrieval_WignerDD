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
addpath ./src/

% For repeatability, set random seed to ...
rng(1234);


%% Parameters
dval =  2^08;                   % signal dimension
rho = round(2.0*log2(dval));  % support of m^hat
addnoise =  true;               % Add noise?
magEst = 'eig';                 % Eigenvector based magnitude estimation
angSync = 'graphleig';          % Graph-Laplacian angulay synchronization


% Print out problem parameters
fprintf( '\n\n --------------------------------------------------- \n' );
fprintf(     '    Phase Retrieval using Wigner Deconvolution \n' );
fprintf(     ' --------------------------------------------------- \n' );

fprintf( ' Problem size, d ~ %d \n', dval );
fprintf( ' Bandwidth of mask (card. of support of m^hat), rho = %d \n', ...
                                     rho );


% Do this for several snr values
snrvals = (10:10:60).';
% Do this for several kappa values
kappavals = (4:4:16);
% ... and no. of trials
ntrials = 100;

% Store error here
error = zeros( length(snrvals), length(kappavals) );

for isnr = 1:length(snrvals)

%% Parameters
snr = snrvals(isnr);            % SNR of noise to be added
    
fprintf( '\n Now running tests at %ddB SNR... \n', snr );
fprintf( ' (Avg.) reconstruction error is \n' );
    
for ikappa = 1:length(kappavals)    
    
% More parameters
kappa = kappavals(ikappa);      % Note: kappa should obey L=rho-1+kappa
L = rho+kappa-1;              % subsampling in space - no. of shifts

% Fix d to satisfy divisibility requirements
d = L*round(dval/L);

fprintf( '   (d,kappa,L) = (%d,%d,%d) ', d,kappa, L );

% Store problem parameters in a structure
pars = struct(  'd', d, 'L', L, 'kappa', kappa, 'pproc', false, ...
                'magEst', magEst, 'angSync', angSync );


%% Choice of mask m

% Note: with random masks, the mask parameter \mu (associated with
% robustness) can be very small (or possibly zero)
% we loop over a few mask choices and choose the best one

% initialization
idx = 1; mu = 0;

while( idx<=100 || mu<=1e-5 )
    % define m^hat
    m_hat = zeros(d,1);
    m_hat(1:rho) = (1+0.5*rand(rho,1)) .* exp(1i*2*pi*rand(rho,1));

    % here is the mask in physical space (for reference)
    m = ifft(m_hat);

    % Pre-computation (of terms involving masks)
    mask_precomp_ld = zeros(d,kappa);
    mask_precomp_ud = zeros(d,kappa-1);

    for alpha = 0:kappa-1
        tmp = fft( m_hat.*circshift(conj(m_hat),-alpha), d );
        mask_precomp_ld(:,alpha+1) = tmp(1:end);
    end

    lb_idx = L-(kappa-1);
    for alpha = L-(kappa-1):L-1
        tmp = fft( m_hat.*circshift(conj(m_hat),L-alpha), d );
        mask_precomp_ud(:,alpha-lb_idx+1) = tmp(1:end);
    end
    
    % mu value
    mu_current = min(abs([mask_precomp_ld(:); mask_precomp_ud(:)]));
    
    % update best choice of masks
    if(mu<mu_current)
        mu = mu_current;
        mhat_best = m_hat;
        m_best = m;
        mask_precomp_ld_best = mask_precomp_ld;
        mask_precomp_ud_best = mask_precomp_ud;
    end
        
    % update loop counter
    idx = idx+1;
end

% Store these in a structure
mask = struct(  'ld', mask_precomp_ld_best, ...
                'ud', mask_precomp_ud_best, ...
                 'm', m_best );

% use the "best" mask
m = m_best;


for itrial = 1:ntrials

%% Choice of signal x
x = randn(d,1)+1i*randn(d,1); % random complex signal x


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


% Solve using Wigner Deconvolution
[xrec, etime] = WignerPR(Y, pars, mask);

% Correct for global phase factor
phaseOffset = angle( (xrec'*x) / (x'*x) );
xrec = xrec * exp(1i*phaseOffset);

% Reconstruction error
errordB = 10*log10( norm(x-xrec)^2 / norm(x)^2 );
error(isnr, ikappa) = error(isnr, ikappa) + errordB/ntrials;

end

fprintf( '%3.2fdB \n', error(isnr,ikappa) );

end

end

% plot
figure; plot(snrvals, error(:,1), '-+', 'linewidth', 2); hold on
plot(snrvals, error(:,2), '--o', 'linewidth', 2);
plot(snrvals, error(:,3), ':x', 'linewidth', 2);
plot(snrvals, error(:,4), '-.s', 'linewidth', 2);
axis([0 65 -65 10]); grid
xlabel( 'Noise Level in SNR (dB)', 'interpreter', 'latex', 'fontsize', 14 )
ylabel( 'Reconstruction Error (in dB)', 'interpreter', ...
                                                'latex', 'fontsize', 14 )
title( 'Robustness to Measurement Noise, $d\approx 256, D=Ld$', ...
                                 'interpreter', 'latex', 'fontsize', 14 )
legend( {'$L=19$', '$L=23$', '$L=27$', '$L=31$'}, ...
                                 'interpreter', 'latex', 'fontsize', 14 )
drawnow;
                             
% Save data
% save( 'noise_L.mat', 'snrvals', 'error', 'kappavals' );
                             