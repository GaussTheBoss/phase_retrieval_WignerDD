% Script to evaluate mask constant mu_2 for different types of masks
%
%% (Ptychographic) Phase Retrieval using Wigner Deconvolution
%
% Script to implement phase retrieval using using Wigner Deconvolution. 
% We solve the phase recovery problem
%
%

clear; close all; clc
addpath ./src/

% For repeatability, set random seed to ...
rng(1234);


% Repeat for several values of delta
deltavals = (2:12).';

% Store mu_1 values here
muvals = zeros(length(deltavals), 4);
ddelta = zeros(length(deltavals), 2);

for idvals = 1:length(deltavals)

    
%% Parameters
d =  2^06;                   % signal dimension (will only choose odd d)

delta = deltavals(idvals);

addnoise =  true;               % Add noise?
magEst = 'eig';                 % Eigenvector based magnitude estimation
angSync = 'graphleig';          % Graph-Laplacian angulay synchronization

% More parameters
kappa = round(delta/1);           % Note: kappa should obey L=delta-1+kappa
K = delta+kappa-1;              % subsampling in space - no. of shifts

% Fix d to satisfy divisibility requirements
d = K*round(d/K);

if(mod(d,2)==0)
    d = d + 1;                      % Only consider odd d
end

% Store these values
ddelta(idvals,1) = d;
ddelta(idvals,2) = delta;

% Store problem parameters in a structure
pars = struct(  'd', d, 'K', K, 'kappa', kappa, 'pproc', false, ...
                'magEst', magEst, 'angSync', angSync );

% Print out problem parameters
fprintf( '\n\n --------------------------------------------------- \n' );
fprintf( ' Problem size, d ~ %d \n', d );
fprintf( ' Support of mask (card. of support of m), delta = %d \n', ...
                                     delta );


%% Choice of mask m - Random mask

% Note: with random masks, the mask parameter \mu (associated with
% robustness) can be very small (or possibly zero)
% we loop over a few mask choices and choose the best one

% initialization
idx = 1; mu = 0;

% Avg mu
mu_avg = 0;

ntrials = 100;


% while( idx<=100 || mu<=1e-5 )
while( idx<=ntrials )

    % define m
    m = zeros(d,1);
    

    % Random masks
    m(1:delta) = (1+0.5*rand(delta,1)) .* exp(1i*2*pi*rand(delta,1));
        
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
    mu_current = min(abs([mask_precomp_ld(:); mask_precomp_ud(:)]));
    
    mu_avg = mu_avg + mu_current/ntrials;
    
    % update best choice of masks
    if(mu<mu_current)
        mu = mu_current;
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

fprintf('\n Average value of mu_2 for random mask is %3.3e\n', mu_avg );





%% Choice of mask m - (Deterministic) Exponential mask

% initialization
idx = 1; mu = 0;

% Avg mu
mu_avg = 0;

ntrials = 1;


% while( idx<=100 || mu<=1e-5 )
while( idx<=ntrials )

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
    mu_current = min(abs([mask_precomp_ld(:); mask_precomp_ud(:)]));
    
    mu_avg = mu_avg + mu_current/ntrials;
    
    % update best choice of masks
    if(mu<mu_current)
        mu = mu_current;
        m_best = m;
        mask_precomp_ld_best = mask_precomp_ld;
        mask_precomp_ud_best = mask_precomp_ud;
    end
        
    % update loop counter
    idx = idx+1;

end


% Display mask-constant mu
mu_avg_rand = mu_avg;
muvals(idvals, 2) = mu_avg;

fprintf(' Average value of mu_2 for (deterministic) exponential mask is %3.3e\n', mu_avg );





%% Choice of mask m - Random Symmetric mask

% Note: with random masks, the mask parameter \mu (associated with
% robustness) can be very small (or possibly zero)
% we loop over a few mask choices and choose the best one

% initialization
idx = 1; mu = 0;

% Avg mu
mu_avg = 0;

ntrials = 100;


% while( idx<=100 || mu<=1e-5 )
while( idx<=ntrials )

    % define m
    m = zeros(d,1);
    

    % Random symmetric (within support) masks
    if( mod(delta,2) == 0 )
        m(1:delta/2) = (1+0.5*rand(delta/2,1)) .* exp(1i*2*pi*rand(delta/2,1));
        m(delta/2+1:delta) = (flipud(m(1:delta/2)));
    else
        m(1:ceil(delta/2)) = (1+0.5*rand(ceil(delta/2),1)) .* exp(1i*2*pi*rand(ceil(delta/2),1));
        m(ceil(delta/2)+1:delta) = (flipud(m(1:floor(delta/2))));
    end

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
    mu_current = min(abs([mask_precomp_ld(:); mask_precomp_ud(:)]));
    
    mu_avg = mu_avg + mu_current/ntrials;
    
    % update best choice of masks
    if(mu<mu_current)
        mu = mu_current;
        m_best = m;
        mask_precomp_ld_best = mask_precomp_ld;
        mask_precomp_ud_best = mask_precomp_ud;
    end
        
    % update loop counter
    idx = idx+1;

end


% Display mask-constant mu
mu_avg_rand = mu_avg;
muvals(idvals, 3) = mu_avg;

fprintf(' Average value of mu_2 for random symmetric mask is %3.3e\n', mu_avg );





%% Choice of mask m - Gaussian mask

% initialization
idx = 1; mu = 0;

% Avg mu
mu_avg = 0;

ntrials = 1;


% while( idx<=100 || mu<=1e-5 )
while( idx<=ntrials )

    % define m
    m = zeros(d,1);
    

    % Gaussian masks
    m(1:delta) = window(@gausswin, delta, 2);
%     m(1:delta) = window(@hamming, delta);
        
    
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
    mu_current = min(abs([mask_precomp_ld(:); mask_precomp_ud(:)]));
    
    mu_avg = mu_avg + mu_current/ntrials;
    
    % update best choice of masks
    if(mu<mu_current)
        mu = mu_current;
        m_best = m;
        mask_precomp_ld_best = mask_precomp_ld;
        mask_precomp_ud_best = mask_precomp_ud;
    end
        
    % update loop counter
    idx = idx+1;

end


% Display mask-constant mu
mu_avg_rand = mu_avg;
muvals(idvals, 4) = mu_avg;

fprintf(' Average value of mu_2 for (deterministic) Gaussian mask is %3.3e\n', mu_avg );


end


                             
