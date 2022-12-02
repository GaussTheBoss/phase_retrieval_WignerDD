% Script to evaluate mask constant mu_1 for different types of masks
%

%% (Ptychographic) Phase Retrieval using Wigner Deconvolution
%
% Script to implement phase retrieval using using Wigner Deconvolution. 
% We solve the phase recovery problem
%

clear; close all; clc
addpath ./src/

% For repeatability, set random seed to ...
rng(1234);


% Repeat for several values of d
rhovals = (2:12).';

% Store mu_1 values here
muvals = zeros(length(rhovals), 4);
ddelta = zeros(length(rhovals), 2);

for idvals = 1:length(rhovals)


%% Parameters
d =  2^06;                      % signal dimension

rho = rhovals(idvals);

kappa = round(rho/1);         % Note: kappa should obey L=delta-1+kappa
L = rho+kappa-1;              % subsampling in space - no. of shifts

% Mask parameters
maskType = 'random';            % random masks


% Noise parameters
addnoise =  true;               % Add noise?

% Postprocessing and other algorithmic parameters
% magEst = 'diag';                % diagonal magnitude estimation
magEst = 'eig';                 % eigenvector based magnitude estimation
angSync = 'graphleig';          % Graph Laplacian angular synchronization 
pproc = false;                  % Post-process using HIO+ER?

% Print mask metrics?
maskmetrics = false;


% Fix d to satisfy divisibility requirements
d = L*round(d/L);


if(mod(d,2)==0)
    d = d + 1;                      % Only consider odd d
end


% Store these values
ddelta(idvals,1) = d;
ddelta(idvals,2) = rho;

% Print out problem parameters
fprintf( '\n\n --------------------------------------------------- \n' );
fprintf( ' Problem size, d ~ %d \n', d );
fprintf( ' Bandwidth of mask (card. of support of m^hat), rho = %d \n', ...
                                     rho );
fprintf( ' No. of shifts, L = %d \n', L );                                 

% Store problem parameters in a structure
pars = struct( 'd', d, 'L', L, 'kappa', kappa, 'pproc', pproc, ...
                        'magEst', magEst, 'angSync', angSync);


%% Choice of mask m

% Note: with random masks, the mask parameter \mu (associated with
% robustness) can be very small (or possibly zero)
% we loop over a few mask choices and choose the best one


%% Random masks

% initialization
idx = 1; mu = 0;

% Avg mu
mu_avg = 0;

ntrials = 100;

% while( idx<=100 || mu<=1e-5 )
while( idx<=ntrials )
    % define m^hat
    m_hat = zeros(d,1);
    
    % Random masks
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
    
    mu_avg = mu_avg + mu_current/ntrials;
    
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

% Display mask-constant mu
mu_avg_rand = mu_avg;
muvals(idvals, 1) = mu_avg;


%% Exponential masks

% initialization
idx = 1; mu = 0;

% Avg mu
mu_avg = 0;

ntrials = 1;

% while( idx<=100 || mu<=1e-5 )
while( idx<=ntrials )
    % define m^hat
    m_hat = zeros(d,1);
    
    % Exp. masks
    a = max(4, (rho-1)/2);
    m_hat(1:rho) = exp(-(0:rho-1).'/a)/((2*rho-1)^.25);
    

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
    
    mu_avg = mu_avg + mu_current/ntrials;
    
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

% Display mask-constant mu
mu_avg_exp = mu_avg;
muvals(idvals, 2) = mu_avg;



%% Random symmetric masks

% initialization
idx = 1; mu = 0;

% Avg mu
mu_avg = 0;

ntrials = 100;

% while( idx<=100 || mu<=1e-5 )
while( idx<=ntrials )
    % define m^hat
    m_hat = zeros(d,1);
    
    % Random symmetric (within support) masks
    if( mod(rho,2) == 0 )
        m_hat(1:rho/2) = (1+0.5*rand(rho/2,1)) .* exp(1i*2*pi*rand(rho/2,1));
        m_hat(rho/2+1:rho) = conj(flipud(m_hat(1:rho/2)));
    else
        m_hat(1:ceil(rho/2)) = (1+0.5*rand(ceil(rho/2),1)) .* exp(1i*2*pi*rand(ceil(rho/2),1));
        m_hat(ceil(rho/2)+1:rho) = conj(flipud(m_hat(1:floor(rho/2))));
    end

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
    
    mu_avg = mu_avg + mu_current/ntrials;
    
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

% Display mask-constant mu
mu_avg_randsym = mu_avg;
muvals(idvals, 3) = mu_avg;



%% Gaussian masks

% initialization
idx = 1; mu = 0;

% Avg mu
mu_avg = 0;

ntrials = 1;

% while( idx<=100 || mu<=1e-5 )
while( idx<=ntrials )
    % define m^hat
    m_hat = zeros(d,1);
    
    % Gaussian masks
    m_hat(1:rho) = window(@gausswin, rho, 2);
%     m_hat(1:delta) = window(@hamming, delta);

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
    
    mu_avg = mu_avg + mu_current/ntrials;
    
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

% Display mask-constant mu
mu_avg_gauss = mu_avg;
muvals(idvals, 4) = mu_avg;


fprintf('\n Average value of mu_1 (for random, exponential, random-symmetric, Gaussina) is %3.3e, %3.3e, %3.3e, %3.3e\n', ...
        mu_avg_rand, mu_avg_exp, mu_avg_randsym, mu_avg_gauss );


               
end
