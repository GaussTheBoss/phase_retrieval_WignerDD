% Copyright (c) 2018- Michigan State University, University of 
% Michigan-Dearborn and the CHARMS Research Group
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

%% (Ptychographic) Phase Retrieval using HIO+ER
%
% Script to implement phase retrieval using using HIO+ER
%   
% For a description of the method, see arXiv preprint
% https://arxiv.org/abs/1907.10773
%

% clear; close all; clc
clear;
addpath ./src/

% For repeatability, set random seed to ...
rng(12345);


%% Parameters
d =  2^08;                      % signal dimension
delta = round(1.25*log2(d));    % support of m^hat
addnoise =  true;               % Add noise?
snr = 30;                       % SNR of noise to be added
kappa = round(delta/1);         % Note: kappa should obey L=delta-1+kappa
L = delta+kappa-1;              % subsampling in space - no. of shifts

% Fix d to satisfy divisibility requirements
d = L*round(d/L);

% some parameter checks
% Require that L|d
if mod(d/L,1)~=0
    disp('Need L|d; Choose a different value for L');
    return
end

% Print out problem parameters
fprintf( '\n\n --------------------------------------------------- \n' );
fprintf(     '    Phase Retrieval using HIO+ER \n' );
fprintf(     ' --------------------------------------------------- \n' );

fprintf( ' Problem size, d = %d \n', d );
fprintf( ' Bandwidth of mask (card. of support of m^hat), delta = %d \n', ...
                                     delta );
fprintf( ' No. of shifts (in space), L = %d \n', L );
fprintf( ' Total no. of measurements, D = %d = %dd\n\n', d*L, L );

if( addnoise )
    fprintf( ' Noisy simulation? - Yes \n' );
    fprintf( ' Added noise (SNR, dB) = %3.2f \n\n', snr );
else
    fprintf( ' Noisy simulation? - No \n\n' );
end


%% Choice of mask m

% define m^hat
m_hat = zeros(d,1);
m_hat(1:delta) = (1+rand(delta,1)) .* exp(1i*2*pi*rand(delta,1));

% here is the mask in physical space (for reference)
m = ifft(m_hat);


% Repeat experiments for 
ntrials = 100;

% Max no of iterations
maxit = 1200;
%(this is done in blocks of 30)
nB = 30;

% Store errors here
error = zeros(maxit/nB, 3);
etime = zeros( size(error) );


for itrials = 1:ntrials
    
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


%% HIO+ER solution

for its = nB:nB:maxit
    tic;
    xrec = randn(d,1);
    xrec = postprocess_HIOER( xrec, Y, m, d, L, [its, 15, 15] );

    % Correct for global phase factor
    phaseOffset = angle( (xrec'*x) / (x'*x) );
    xrec = xrec * exp(1i*phaseOffset);

    % Store error and execution times
    etime(its/nB,1) = etime(its/nB,1) + toc/ntrials;
    
    errordB = 10*log10( norm(x-xrec)^2 / norm(x)^2 );
    error(its/nB,1) = error(its/nB,1) + errordB/ntrials;
    
    
    
    % Now do this for another combination of (HIO,ER) iterations
    tic;
    xrec = randn(d,1);
    xrec = postprocess_HIOER( xrec, Y, m, d, L, [its, 20, 10] );

    % Correct for global phase factor
    phaseOffset = angle( (xrec'*x) / (x'*x) );
    xrec = xrec * exp(1i*phaseOffset);

    % Store error and execution times
    etime(its/nB,2) = etime(its/nB,2) + toc/ntrials;
    
    errordB = 10*log10( norm(x-xrec)^2 / norm(x)^2 );
    error(its/nB,2) = error(its/nB,2) + errordB/ntrials;
    
    
    
    % Now do this for another combination of (HIO,ER) iterations
    tic;
    xrec = randn(d,1);
    xrec = postprocess_HIOER( xrec, Y, m, d, L, [its, 25, 5] );

    % Correct for global phase factor
    phaseOffset = angle( (xrec'*x) / (x'*x) );
    xrec = xrec * exp(1i*phaseOffset);

    % Store error and execution times
    etime(its/nB,3) = etime(its/nB,3) + toc/ntrials;
    
    errordB = 10*log10( norm(x-xrec)^2 / norm(x)^2 );
    error(its/nB,3) = error(its/nB,3) + errordB/ntrials;
    
end

end

%% Plot error vs iteration count

figure; plot( (nB:nB:maxit).', error(:,1), '-o', 'linewidth', 2 ); hold on
plot( (nB:nB:maxit).', error(:,2), '-+', 'linewidth', 2 );
plot( (nB:nB:maxit).', error(:,3), '-s', 'linewidth', 2 );
xlabel( 'No. of iterations', 'interpreter', 'latex', 'fontsize', 14 ); 
ylabel( 'Reconstruction Error (dB)', 'interpreter', 'latex', 'fontsize', 14 );
title( 'Error vs Iteration Count (HIO+ER, 30dB noise)', ...
                            'interpreter', 'latex', 'fontsize', 14 ); 
legend( {'(HIO,ER) = (15,15)', '(HIO,ER) = (20,10)', '(HIO,ER) = (25,5)'}, ...
                        'interpreter', 'latex', 'fontsize', 14 ); 
grid; drawnow;

% Save data
% save( 'HIO_error_vs_iter.mat', 'error', 'etime', 'nB', 'maxit' );