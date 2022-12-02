% Copyright (c) 2018- Michigan State University, 
% University of Michigan-Dearborn and the CHARMS Research Group
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


% Run execution time tests for 
dvals = 2.^(4:10).';
nPL = 135;               % cutoff for PhaseLift
 
% ... and no. of trials
ntrials = 100;

% Store error and execution times here
error = zeros( length(dvals), 5 );
etime = zeros( size(error) );

for id = 1:length(dvals)

%% Parameters
dv =  dvals(id);                % signal dimension
delta = round(1.25*log2(dv));   % support of m^hat
kappa = round(delta/2);         % Note: kappa should obey L=delta-1+kappa
K = delta+kappa-1;              % subsampling in space - no. of shifts

% Fix d to satisfy divisibility requirements
d = K*round(dv/K);
dvals(id) = d;

addnoise =  true;               % Add noise?
snr = 20;                       % SNR of noise to be added
pproc =  false;                 % Post-process using HIO+ER?

fprintf( '\n\n Now running tests at d=%d... (delta,kappa,K)=(%d,%d,%d)', ...
                                    d, delta, kappa, K );

% Store problem parameters in a structure
pars = struct( 'd', d, 'K', K, 'kappa', kappa, ...
                'pproc', pproc, 'snr', snr );


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


if(d<=nPL)
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

end


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
% Y_alt = spectrogram(x, m, rho-1, d);

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
    

%% Solve using Wigner Deconvolution (without post-processing)
pars.pproc = false; pars.magEst = 'diag'; pars.angSync = 'topeig';
[xrec, et] = WignerPR(Y, pars, mask);

% Correct for global phase factor
phaseOffset = angle( (xrec'*x) / (x'*x) );
xrec = xrec * exp(1i*phaseOffset);

% Reconstruction error
errordB = 10*log10( norm(x-xrec)^2 / norm(x)^2 );
error(id, 1) = error(id, 1) + errordB/ntrials;
etime(id, 1) = etime(id, 1) + et/ntrials;


%% Solve using Wigner Deconvolution (with post-processing)
pars.pproc = false; pars.magEst = 'eig'; pars.angSync = 'graphleig';
[xrec, et] = WignerPR(Y, pars, mask);

% Correct for global phase factor
phaseOffset = angle( (xrec'*x) / (x'*x) );
xrec = xrec * exp(1i*phaseOffset);

% Reconstruction error
errordB = 10*log10( norm(x-xrec)^2 / norm(x)^2 );
error(id, 2) = error(id, 2) + errordB/ntrials;
etime(id, 2) = etime(id, 2) + et/ntrials;


%% PhaseLift solution
% Not sure why this works as a stand-alone script but not here
if(d<=nPL)
[xrec, et] = PhaseLiftPR(Y(:), pars, tfocsLinOp);

% Correct for global phase factor
phaseOffset = angle( (xrec'*x) / (x'*x) );
xrec = xrec * exp(1i*phaseOffset);

% Reconstruction error
errordB = 10*log10( norm(x-xrec)^2 / norm(x)^2 );
error(id, 3) = error(id, 3) + errordB/ntrials;
etime(id, 3) = etime(id, 3) + et/ntrials;
end

%% Add WirtingerFlow implementation here
[xrec, et] = WF_EigPR( delta, d, x, snr, pars_WF );

% Correct for global phase factor
phaseOffset = angle( (xrec'*x) / (x'*x) );
xrec = xrec * exp(1i*phaseOffset);

% Reconstruction error
errordB = 10*log10( norm(x-xrec)^2 / norm(x)^2 );
error(id, 4) = error(id, 4) + errordB/ntrials;
etime(id, 4) = etime(id, 4) + et/ntrials;


%% Add HIO+ER implementation here
xrec = zeros(d,1);
tic; 
xrec = HIO_BlockPR_EigPR( xrec, measurements, maskEntries, pars_HIOER );
et = toc;

% Correct for global phase factor
phaseOffset = angle( (xrec'*x) / (x'*x) );
xrec = xrec * exp(1i*phaseOffset);                                   

% Reconstruction error
errordB = 10*log10( norm(x-xrec)^2 / norm(x)^2 );
error(id, 5) = error(id, 5) + errordB/ntrials;
etime(id, 5) = etime(id, 5) + et/ntrials;


end

fprintf( '\n  Exec. time using Alg. 1 is %3.2fs', ...
                                    etime(id, 1) );
fprintf( '\n  Exec. time using Alg. 1 (w/ improved mag. est., ang. sync.) is %3.2fs', ...
                                    etime(id, 2) );
if(d<=nPL)
    fprintf( '\n  Exec. time using PhaseLift is %3.2fs', etime(id, 3) );
end                                
fprintf( '\n  Exec. time using WirtingerFlow is %3.2fs', etime(id, 4) );
fprintf( '\n  Exec. time using HIO+ER is %3.2fs', etime(id, 5) );
                                

end


% plot
% Note: you may need to change the axis limits depending on the speed of 
% your hardware/software implementation
figure; 
loglog(dvals(dvals<=nPL), etime(dvals<=nPL,3), '-x', 'linewidth', 2); hold on
loglog(dvals, etime(:,4), ':d', 'linewidth', 2);
loglog(dvals, etime(:,5), '-.s', 'linewidth', 2);
loglog(dvals, etime(:,2), '-+', 'linewidth', 2);
loglog(dvals, etime(:,1), '--o', 'linewidth', 2);
xlim([8 2*2^10]); ylim([3e-5 1e3]); grid
ylabel( 'Execution Time (in secs.)', 'interpreter', 'latex', 'fontsize', 14 )
xlabel( 'Problem Size, $d$', 'interpreter', 'latex', 'fontsize', 14 )
title( 'Execution Time', 'interpreter', 'latex', 'fontsize', 14 )
legend( {'PhaseLift', 'WirtingerFlow', 'HIO+ER', ...
            'Lemma 11 (w/ Eig. Mag. Est., Graph Ang. Sync.)', ...
            'Lemma 11 Based Approach' }, ...
            'location', 'southeast', 'interpreter', 'latex', 'fontsize', 14 )


dt = linspace(3e1,9e1,50);                             
loglog( dt, 3e-5*(dt.^3), 'k--', 'linewidth', 2 ,'HandleVisibility','off'); 
dt = linspace(3e1,2e2,50); 
loglog( dt, 2e-3*dt.*log2(dt), 'k--', 'linewidth', 2 ...
                                    ,'HandleVisibility','off');
drawnow;
                                
% Save data
% save( 'etime_lemma11.mat', 'dvals', 'etime', 'error' );

fprintf( '\n\n ' );
