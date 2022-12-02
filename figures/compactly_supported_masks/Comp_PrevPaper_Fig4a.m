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
rng(12345);


%% Parameters
d =  2^08;                      % signal dimension
rho = round(1.25*log2(d));      % support of m
addnoise =  true;               % Add noise?
snr = 50;                       % SNR of noise to be added
kappa = round(rho/1);           % Note: kappa should obey L=rho-1+kappa
K = rho+kappa-1;                % subsampling in frequencies - no. of modes

angSync = 'topeig';             % Use top eig. of normalized matrix
magEst = 'diag';                % Use sqrt. of diagonals for magnitudes
pproc =  false;                 % Post-process using HIO+ER?

% Fix d to satisfy divisibility requirements
d = K*round(d/K);

% some parameter checks
% Require that K|d
if mod(d/K,1)~=0
    disp('Need K|d; Choose a different value for K');
    return
end

% Print out problem parameters
fprintf( '\n\n --------------------------------------------------- \n' );
fprintf(     '    Phase Retrieval using Wigner Deconvolution \n' );
fprintf(     ' --------------------------------------------------- \n' );

fprintf( ' Problem size, d = %d \n', d );
fprintf( ' Support parameter (card. of support of m), rho = %d \n', rho );
fprintf( ' No. of modes (frequencies), K = %d \n', K );
fprintf( ' Total no. of measurements, D = %d = %dd\n\n', d*K, K );

if( addnoise )
    fprintf( ' Noisy simulation? - Yes \n' );
    fprintf( ' Added noise (SNR, dB) = %3.2f \n\n', snr );
else
    fprintf( ' Noisy simulation? - No \n\n' );
end

% Store problem parameters in a structure
pars = struct( 'd', d, 'K', K, 'kappa', kappa, ...
               'angSync', angSync, 'magEst', magEst, 'pproc', pproc);


           
%% These are for the implementation from [33] (EigPR)
% HIO+ER iterations
niter = 600;
npp_iter = 60;
nHIO = 20;
nER = 10;


% Store the above parameters in a structure
    pars_EigPR_Alg1 = struct(      ...
                        'd',                d, ...
                        ...
                        'addNoise',         addnoise, ...
                        'snr',              snr, ...
                        ...
                        'testFncType',      'randGauss', ...
                        'useFlattening',    false, ...
                        ...
                        'maskType',         'fourier', ...
                        'delta',            rho, ...
                        'over',             1, ...
                        'fullMat',          false, ...                        
                        ...
                        'alg',              'BlockPR', ...
                        'angSyncMethod',    'eig', ...
                        'MagEstType',       'diag', ...                        
                        'altProj',          false, ...
                        'altProjIters',     100, ...
                        'threshold',        1e-10 ...
                        ...
                        );

           
% BlockPR - from [33], Alg. 1 w/ Mag. Est + Alt. Proj (ER)
    pars_EigPR_PP = struct(      ...
                        'd',                d, ...
                        ...
                        'addNoise',         addnoise, ...
                        'snr',              snr, ...
                        ...
                        'testFncType',      'randGauss', ...
                        'useFlattening',    false, ...
                        ...
                        'maskType',         'fourier', ...
                        'delta',            rho, ...
                        'over',             1, ...
                        'fullMat',          false, ...                        
                        ...
                        'alg',              'BlockPR', ...
                        'angSyncMethod',    'eig', ...
                        'MagEstType',       'eig', ...
                        'altProj',          true, ...                        
                        'altProjIters',     npp_iter, ...
                        'niter_HIO',        nHIO, ...
                        'niter_ER',         nER, ...                         
                        'threshold',        1e-10 ...
                        ...
                        );                    
                    
           
%% Choice of mask m

% define m
m = zeros(d,1);

% Exponential mask (deterministic, real)
% parameter
a = max(4, (rho-1)/2);
m(1:rho) = exp(-(0:rho-1).'/a)/((2*rho-1)^.25);


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

% Store these in a structure
mask = struct( 'ld', mask_precomp_ld, 'ud', mask_precomp_ud, 'm', m );


% Do this for several snr values
snrvals = (10:10:60).';
% ... and no. of trials
ntrials = 100;

% Store error here
error = zeros( length(snrvals), 4 );

for isnr = 1:length(snrvals)

%% Parameters
snr = snrvals(isnr);            % SNR of noise to be added
    
fprintf( '\n Now running tests at %ddB SNR... \n', snr );
    

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


%% For the implementations from [33] (EigPR)

% Construct masks for generating measurements
[maskEntries, blockDiag] = constructMask_EigPR(pars_EigPR_Alg1);

% Generate measurements
% TODO: should be able to use the measurements from above
pars_EigPR_Alg1.snr = snr;
pars_EigPR_PP.snr = snr;
[measurements, pars_EigPR_Alg1] = ...
        generateMeasurements_EigPR( x, maskEntries, pars_EigPR_Alg1 );
% Record noise power for other algorithms
pars_EigPR_PP.noise_power = pars_EigPR_Alg1.noise_power;



%% Solve using Wigner Deconvolution (no post-processing)
pars.magEst = 'diag';
pars.angSync = 'topeig';
[xrec, etime] = WignerPR(Y, pars, mask);

% Correct for global phase factor
phaseOffset = angle( (xrec'*x) / (x'*x) );
xrec = xrec * exp(1i*phaseOffset);

% Reconstruction error
errordB = 10*log10( norm(x-xrec)^2 / norm(x)^2 );
error(isnr, 1) = error(isnr, 1) + errordB/ntrials;


%% Solve using Wigner Deconvolution (with post-processing)
pars.magEst = 'eig';
pars.angSync = 'graphleig';
[xrec, etime] = WignerPR(Y, pars, mask);

% Correct for global phase factor
phaseOffset = angle( (xrec'*x) / (x'*x) );
xrec = xrec * exp(1i*phaseOffset);

% Reconstruction error
errordB = 10*log10( norm(x-xrec)^2 / norm(x)^2 );
error(isnr, 2) = error(isnr, 2) + errordB/ntrials;



%% Solve using Alg. 1 from [33] (EigPR, no post-processing)
[xrec, etime] = blockPR_EigPR( measurements, maskEntries, blockDiag, pars_EigPR_Alg1);

% Correct for global phase factor
phaseOffset = angle( (xrec'*x) / (x'*x) );
xrec = xrec * exp(1i*phaseOffset);

% Reconstruction error
errordB = 10*log10( norm(x-xrec)^2 / norm(x)^2 );
error(isnr, 3) = error(isnr, 3) + errordB/ntrials;


%% Solve using Alg. 1 from [33] (EigPR, w/ post-processing)
[xrec, etime] = blockPR_EigPR( measurements, maskEntries, blockDiag, pars_EigPR_PP);

% Correct for global phase factor
phaseOffset = angle( (xrec'*x) / (x'*x) );
xrec = xrec * exp(1i*phaseOffset);

% Reconstruction error
errordB = 10*log10( norm(x-xrec)^2 / norm(x)^2 );
error(isnr, 4) = error(isnr, 4) + errordB/ntrials;


end

fprintf( ' Reconstruction error (Lemma 11) is %3.2fdB\n', error(isnr,1) );
fprintf( ' Reconstruction error w/ Eig. Mag. Est., Graph L. Ang. Sync. is %3.2fdB\n', error(isnr,2) );
fprintf( ' Reconstruction error from [33], Alg. 1 is %3.2fdB\n', error(isnr,3) );
fprintf( ' Reconstruction error from [33], Alg. 1 w/ HIO+ER post-proc., Eig. Mag. Est. is %3.2fdB\n', error(isnr,4) );


end

% plot
figure; plot(snrvals, error(:,1), '-+', 'linewidth', 2); hold on
plot(snrvals, error(:,3), '--o', 'linewidth', 2);
plot(snrvals, error(:,4), '--x', 'linewidth', 2);
plot(snrvals, error(:,2), '-s', 'linewidth', 2);
axis([0 65 -70 10]); grid
xlabel( 'Noise Level in SNR (dB)', 'interpreter', 'latex', 'fontsize', 14 )
ylabel( 'Reconstruction Error (in dB)', 'interpreter', ...
                                                'latex', 'fontsize', 14 )
title( 'Robustness to Measurement Noise, $d=247, D=19d$', ...
                                 'interpreter', 'latex', 'fontsize', 14 )
legend( {'Lemma 11', 'from [33], Alg. 1', 'from [33], Alg. 1 w/ post proc.', ...
         'Lemma 11 w/ Eig. Mag. Est., Graph Ang. Sync.'}, ...
     'interpreter', 'latex', 'fontsize', 14, 'location', 'southwest' )
drawnow;
 
% % Save data
% save( 'noise_postproc_K.mat', 'snrvals', 'error' );
                                 