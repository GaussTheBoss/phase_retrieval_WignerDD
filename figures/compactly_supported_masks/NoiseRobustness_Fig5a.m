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
% Change this path if necessary
addpath ../third/TFOCS-1.4.1

% For repeatability, set random seed to ...
rng(1234);


%% Parameters
dval =  2^06;                   % signal dimension
delta = round(1.25*log2(dval)); % support of m^hat
addnoise =  true;               % Add noise?
magEst = 'eig';                 % Eigenvector based magnitude estimation
angSync = 'graphleig';          % Graph-Laplacian angulay synchronization

% More parameters
kappa = round(delta/1);           % Note: kappa should obey L=delta-1+kappa
K = delta+kappa-1;              % subsampling in space - no. of shifts

% Fix d to satisfy divisibility requirements
d = K*round(dval/K);

% Store problem parameters in a structure
pars = struct(  'd', d, 'K', K, 'kappa', kappa, 'pproc', false, ...
                'magEst', magEst, 'angSync', angSync );

% Print out problem parameters
fprintf( '\n\n --------------------------------------------------- \n' );
fprintf(     '    Phase Retrieval using Wigner Deconvolution \n' );
fprintf(     ' --------------------------------------------------- \n' );

fprintf( ' Problem size, d ~ %d \n', dval );
fprintf( ' Support of mask (card. of support of m), delta = %d \n', ...
                                     delta );


% Do this for several snr values
snrvals = (10:10:60).';

% ... and no. of trials
ntrials = 100;

% Store error here
error = zeros( length(snrvals), 5 );

for isnr = 1:length(snrvals)

%% Parameters
snr = snrvals(isnr);            % SNR of noise to be added
pars.snr = snr;
    
fprintf( '\n Now running tests at %ddB SNR... \n', snr );
    

%% Choice of mask m

% use exponential masks from previous paper
% define m
m = zeros(d,1);

% Exponential mask (deterministic, real)
% parameter
a = max(4, (delta-1)/2);
m(1:delta) = exp(-(0:delta-1).'/a)/((2*delta-1)^.25);

% Pre-computation (of terms involving masks)
mask_precomp_ud = zeros(d,kappa);
mask_precomp_ld = zeros(d,kappa-1);

for omega = 0:kappa-1
    mask_precomp_ud(:,omega+1) = fft( m.*circshift(conj(m),-omega), d );
end

lb_idx = K-(kappa-1);
for omega = K-(kappa-1):K-1
    mask_precomp_ld(:,omega-lb_idx+1) = fft( m.*circshift(conj(m),K-omega), d );
end

mask_precomp_ld = [mask_precomp_ld(1,:); flipud(mask_precomp_ld(2:end,:))];
mask_precomp_ud = [mask_precomp_ud(1,:); flipud(mask_precomp_ud(2:end,:))];
    
% mu value
mu = min(abs([mask_precomp_ld(:); mask_precomp_ud(:)]));
    

% Store these in a structure
mask = struct(  'ld', mask_precomp_ld, ...
                'ud', mask_precomp_ud, ...
                 'm', m );


% TFOCS requires special initialization of the measurement model
% Here is the measurement matrix
measurementMat = zeros(K*d,d);

% subsampled Fourier matrix
Fd = dftmtx(d);
F_K = Fd(1:d/K:end, :);
clear Fd;

for l=0:d-1
    measurementMat( l*K+1:(l+1)*K, :) = F_K*diag(circshift(m,l));
end

tfocsLinOp = initializeLinopPR( {[d, d], [K*d, 1]}, measurementMat ); 

clear measurementMat;


% For fast implementation of HIO+ER
% Store the above parameters in a structure
    pars_HIOER = struct(      ...
                        'd',                d, ...
                        ...
                        'addNoise',         addnoise, ...
                        'snr',              snr, ...
                        ...
                        'testFncType',      'randGauss', ...
                        'useFlattening',    false, ...
                        ...
                        'maskType',         'fourier', ...
                        'delta',            delta, ...
                        'over',             1, ...
                        'fullMat',          false, ...                        
                        ...
                        'alg',              'BlockPR', ...
                        'angSyncMethod',    'eig', ...
                        'MagEstType',       'diag', ...                        
                        'altProj',          true, ...
                        'threshold',        1e-10, ...
                        'altProjIters',     600, ...
                        'niter_HIO',        25, ...
                        'niter_ER',         5 ...                         
                        ...
                        );

% Construct masks for generating measurements
[maskEntries, blockDiag] = constructMask_EigPR(pars_HIOER);

% For WirtingerFlow
pars_WF = pars;
pars_WF.maskType = 'local';
pars_WF.delta = delta; 
pars_WF.over = 1;
pars_WF.alg = 'EigPR';


for itrial = 1:ntrials

%% Choice of signal x
x = randn(d,1)+1i*randn(d,1); % random complex signal x


%% Measurements
Y = zeros(K,d);             % matrix of (spectrogram) measurements

% Each column   of Y corresponds of a physical shift \ell
% Each row      of Y corresponds to a Fourier  mode  \omega
% TODO: <change this to only compute subsampled measurements>
for l=0:d-1
    tmp = abs( fft(x.*circshift(m,l), d) ).^2;
    Y(:,l+1) = tmp(1:d/K:end);
end

% TODO:
% do this using spectrogram

% noise
if(addnoise)
    signal_power = norm( Y(:) )^2 / (d*K);
    noise_power = signal_power / ( 10^(snr/10) );    
    % Add (real) Gaussian noise of desired variance
    noise = sqrt(noise_power)*randn( size(Y) );
    Y = Y + noise;    
else
    noise_power = 0;
end


% For the HIO+ER implementation
% TODO: should be able to use the measurements above
[measurements, pars_HIOER] = ...
        generateMeasurements_EigPR( x, maskEntries, pars_HIOER );


%% Solve using Wigner Deconvolution (w/ improved mag. est., ang. sync.)
pars.magEst = 'eig';
pars.angSync = 'graphleig'; 
[xrec, etime] = WignerPR(Y, pars, mask);

% Correct for global phase factor
phaseOffset = angle( (xrec'*x) / (x'*x) );
xrec = xrec * exp(1i*phaseOffset);

% Reconstruction error
errordB = 10*log10( norm(x-xrec)^2 / norm(x)^2 );
error(isnr, 1) = error(isnr, 1) + errordB/ntrials;



%% Solve using Wigner Deconvolution (Alg. 1)
pars.magEst = 'diag';
pars.angSync = 'topeig'; 
[xrec, etime] = WignerPR(Y, pars, mask);

% Correct for global phase factor
phaseOffset = angle( (xrec'*x) / (x'*x) );
xrec = xrec * exp(1i*phaseOffset);

% Reconstruction error
errordB = 10*log10( norm(x-xrec)^2 / norm(x)^2 );
error(isnr, 2) = error(isnr, 2) + errordB/ntrials;



%% Solve using PhaseLift
[xrec, etime] = PhaseLiftPR(Y(:), pars, tfocsLinOp);

% Correct for global phase factor
phaseOffset = angle( (xrec'*x) / (x'*x) );
xrec = xrec * exp(1i*phaseOffset);

% Reconstruction error
errordB = 10*log10( norm(x-xrec)^2 / norm(x)^2 );
error(isnr, 3) = error(isnr, 3) + errordB/ntrials;



%% Solve using HIO+ER
[xrec] = HIO_BlockPR_EigPR( zeros(d,1), measurements, maskEntries, pars_HIOER );

% Correct for global phase factor
phaseOffset = angle( (xrec'*x) / (x'*x) );
xrec = xrec * exp(1i*phaseOffset);

% Reconstruction error
errordB = 10*log10( norm(x-xrec)^2 / norm(x)^2 );
error(isnr, 4) = error(isnr, 4) + errordB/ntrials;



%% Solve using WirtingerFlow
[xrec, etime] = WF_EigPR( delta, d, x, snr, pars_WF );

% Correct for global phase factor
phaseOffset = angle( (xrec'*x) / (x'*x) );
xrec = xrec * exp(1i*phaseOffset);

% Reconstruction error
errordB = 10*log10( norm(x-xrec)^2 / norm(x)^2 );
error(isnr, 5) = error(isnr, 5) + errordB/ntrials;



end

fprintf( ' (Avg.) reconstruction error using Alg. 1 is %3.3f\n', error(isnr, 2) );
fprintf( ' (Avg.) reconstruction error using Alg. 1 (w/ improved mag. est., ang. sync. is %3.3f\n', error(isnr, 1) );
fprintf( ' (Avg.) reconstruction error using PhaseLift is %3.3f\n', error(isnr, 3) );
fprintf( ' (Avg.) reconstruction error using HIO+ER is %3.3f\n', error(isnr, 4) );
fprintf( ' (Avg.) reconstruction error using WirtingerFlow is %3.3f\n', error(isnr, 5) );


end


% plot
figure; plot(snrvals, error(:,1), '-+', 'linewidth', 2); hold on
plot(snrvals, error(:,2), '--o', 'linewidth', 2);
plot(snrvals, error(:,3), '-x', 'linewidth', 2);
plot(snrvals, error(:,4), '-.s', 'linewidth', 2);
plot(snrvals, error(:,5), ':d', 'linewidth', 2);
axis([0 65 -80 10]); grid
xlabel( 'Noise Level in SNR (dB)', 'interpreter', 'latex', 'fontsize', 14 )
ylabel( 'Reconstruction Error (in dB)', 'interpreter', ...
                                                'latex', 'fontsize', 14 )
title( 'Robustness to Measurement Noise, $d=60, \delta=8, K=2\delta-1$', ...
                                 'interpreter', 'latex', 'fontsize', 14 )
legend( {'Lemma 11 Based Approach (w/ Eig. Mag. Est., Graph Ang. Sync.', ...
         'Lemma 11 Based Approach', ...
         'PhaseLift', 'HIO+ER', 'WirtingerFlow'}, ...
         'location', 'southwest', 'interpreter', 'latex', 'fontsize', 14 )
drawnow;
     
% Save data
% save( 'noise_lemma11.mat', 'snrvals', 'error' );
                             