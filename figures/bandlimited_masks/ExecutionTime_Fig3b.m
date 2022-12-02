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
rho = round(1.25*log2(dv));   % support of m^hat
kappa = round(rho/2);         % Note: kappa should obey L=rho-1+kappa
L = rho+kappa-1;              % subsampling in space - no. of shifts

% Fix d to satisfy divisibility requirements
d = L*round(dv/L);
dvals(id) = d;

addnoise =  true;               % Add noise?
snr = 40;                       % SNR of noise to be added
pproc =  false;                 % Post-process using HIO+ER?

fprintf( '\n\n Now running tests at d=%d... (rho,kappa,L)=(%d,%d,%d)', ...
                                    d, rho, kappa, L );

% Store problem parameters in a structure
pars = struct( 'd', d, 'L', L, 'kappa', kappa, ...
                'pproc', pproc, 'snr', snr );


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


if(d<=nPL)
% Measurement matrix for PhaseLift implementation
% note - this can be made more efficient, but not with current
% implementation of CVX
measurementMat = zeros(L*d,d);

for l=0:L-1
    shft = l*d/L;
    measurementMat( l*d+1:(l+1)*d ,:) = dftmtx(d)*diag(circshift(m,shft));
end
end

    
for itrial = 1:ntrials

%% Choice of signal x
x = randn(d,1)+1i*randn(d,1); % random complex signal x


%% Measurements
Y = zeros(d,L);             % matrix of (spectrogram) measurements

% Each column   of Y corresponds of a physical shift \ell
% Each row      of Y corresponds to a Fourier  mode  \omega
% TODO: <change this to only compute subsampled measurements>
for l=0:d/L:d-1
    Y(:,l/(d/L)+1) = abs( fft(x.*circshift(m,l), d) ).^2;
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
[xrec, et] = PhaseLiftPR(Y(:), pars, measurementMat);

% Correct for global phase factor
phaseOffset = angle( (xrec'*x) / (x'*x) );
xrec = xrec * exp(1i*phaseOffset);

% Reconstruction error
errordB = 10*log10( norm(x-xrec)^2 / norm(x)^2 );
error(id, 3) = error(id, 3) + errordB/ntrials;
etime(id, 3) = etime(id, 3) + et/ntrials;
end

%% Add WirtingerFlow implementation here
[xrec, et] = WFlowPR(Y, pars, mask);

% Correct for global phase factor
phaseOffset = angle( (xrec'*x) / (x'*x) );
xrec = xrec * exp(1i*phaseOffset);

% Reconstruction error
errordB = 10*log10( norm(x-xrec)^2 / norm(x)^2 );
error(id, 4) = error(id, 4) + errordB/ntrials;
etime(id, 4) = etime(id, 4) + et/ntrials;


%% Add HIO+ER implementation here
xrec = zeros(d,1);
[xrec, et] = postprocess_HIOER( xrec, Y, m, d, L, [400, 25, 5] );

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
loglog(dvals(dvals<=nPL), etime(dvals<=nPL,3), ':o', 'linewidth', 2); hold on
loglog(dvals, etime(:,4), '-.x', 'linewidth', 2);
loglog(dvals, etime(:,5), '--s', 'linewidth', 2);
loglog(dvals, etime(:,2), '-+', 'linewidth', 2);
loglog(dvals, etime(:,1), '-d', 'linewidth', 2);
xlim([8 2*2^10]); ylim([3e-5 5e2]); grid
ylabel( 'Execution Time (in secs.)', 'interpreter', 'latex', 'fontsize', 14 )
xlabel( 'Problem Size, $d$', 'interpreter', 'latex', 'fontsize', 14 )
title( 'Execution Time', 'interpreter', 'latex', 'fontsize', 14 )
legend( {'PhaseLift', 'WirtingerFlow', 'HIO+ER', ...
            'Alg. 1 (w/ modified Steps (5), (7))', 'Alg. 1' }, ...
            'location', 'southeast', 'interpreter', 'latex', 'fontsize', 14 )


dt = linspace(3e1,9e1,50);                             
loglog( dt, 3e-5*(dt.^3), 'k--', 'linewidth', 2 ,'HandleVisibility','off'); 
dt = linspace(3e1,2e2,50); 
loglog( dt, 8e-4*dt.*log2(dt), 'k--', 'linewidth', 2 ...
                                    ,'HandleVisibility','off');
drawnow;
                                
% Save data
% save( 'etime.mat', 'dvals', 'etime', 'error' );

fprintf( '\n\n ' );
