function [ recoveredSoln, etime ] = WF_EigPR( delta, d, trueSig, snr, pars )

% Implements the Wirtinger Flow phase retrieval algorithm
% this function largely uses code by 
%   Mahdi Soltanolkotabi (can be found at)
%       http://www-bcf.usc.edu/~soltanol/WFcode.html
%   Yuxin Chen
%       http://web.stanford.edu/~yxchen/TWF/code.html
%

%% Initialization
% Construct masks for generating measurements

if(strcmp(pars.maskType, 'CDP'))
    % CDP measurements
    L = 2*delta-1;                   % Number of masks  

    % Sample phases: each symbol in alphabet {1, -1, i , -i} has equal prob. 
    Masks = randsrc(d,L,[1i -1i 1 -1]);

    % Sample magnitudes and make masks 
    temp = rand(size(Masks));
    Masks = Masks .* ( (temp <= 0.2)*sqrt(3) + (temp > 0.2)/sqrt(2) );

    % Make linear operators; A is forward map and At its scaled adjoint (At(Y)*numel(Y) is the adjoint) 
    A_CDP = @(I)  fft(conj(Masks) .* repmat(I,[1 L]));  % Input is n x 1 signal, output is n x L array
    At_CDP = @(Y) mean(Masks .* ifft(Y), 2);            % Input is n x L array, output is n x 1 signal          

    
    % measurements
    Y_CDP = abs(A_CDP(trueSig)).^2;  

    % Calculate signal power
    % This is just the norm^2 of the measurements divided by dimension
    signal_power = norm( Y_CDP(:) )^2 / numel(Y_CDP);

    % Calculate noise power to be added
    noise_power = signal_power / ( 10^(snr/10) );
    noise = sqrt(noise_power)*randn( size(Y_CDP) );
    Y_CDP = Y_CDP + noise;


    %% Wirtinger Flow -- CDP 

    tic;

    % Initialization
    npower_iter = 50;                          % Number of power iterations 
    z0 = randn(d,1); z0 = z0/norm(z0,'fro'); % Initial guess 
    for tt = 1:npower_iter, 
        z0 = At_CDP(Y_CDP.*A_CDP(z0)); z0 = z0/norm(z0,'fro');
    end

    normest = sqrt(sum(Y_CDP(:))/numel(Y_CDP)); % Estimate norm to scale eigenvector  
    z = normest * z0;                   % Apply scaling 


    %% Loop

    T = 2500;                           % Max number of iterations
    tau0 = 330;                         % Time constant for step size
    mu = @(t) min(1-exp(-t/tau0), 0.2); % Schedule for step size

    for t = 1:T,
        Bz = A_CDP(z);
        C  = (abs(Bz).^2-Y_CDP) .* Bz;
        grad = At_CDP(C);                    % Wirtinger gradient
        z = z - mu(t)/normest^2 * grad;  % Gradient update 
    end

    % Output quantities
    recoveredSoln = z;
    etime = toc;


elseif(strcmp(pars.maskType, 'local'))

    % Initialization
    
    pars_mask = pars;           % temporary var
    pars_mask.maskType = 'rayan';
    pars_mask.fullMat = true;
    
    % measurement matrix
    [~, ~, measurementMat] = constructMask_EigPR(pars_mask);
    
    A = sparse(measurementMat);
    % measurements
    y = abs(A*trueSig).^2;  

    % Calculate signal power
    % This is just the norm^2 of the measurements divided by dimension
    signal_power = norm( y(:) )^2 / numel(y);

    % Calculate noise power to be added
    noise_power = signal_power / ( 10^(snr/10) );
    noise = sqrt(noise_power)*randn( size(y) );
    y = y + noise;
    m = length(y);

    tic; 
    
    % Spectral initialization
    npower_iter = 1000;                         % Number of power iterations 
    z0 = randn(d,1); z0 = z0/norm(z0,'fro');    % Initial guess 
    for tt = 1:npower_iter                      % Power iterations
        z0 = A'*(y.* (A*z0)); z0 = z0/norm(z0,'fro');
    end

    normest = sqrt(sum(y)/numel(y));    % Estimate norm to scale eigenvector  
    z = normest * z0;                   % Apply scaling 


    %% Loop

    % T = 2500;                           % Max number of iterations
%     T = 7500;                           % Max number of iterations
    
    switch snr
        case 10
            T = 1750;
        case 20
            T = 3250;
        case 30
            T = 4250;
        case 40
            T = 6000;
        case 50
            T = 7500;
        case 60
            T = 10000;
    end
    
    
    % tau0 = 330;                         % Time constant for step size
    tau0 = 1;

    % mu = @(t) min(1-exp(-t/tau0), 0.2); % Schedule for step size
    mu = @(t) min(1-exp(-t/tau0), 10*pars_mask.delta); % Schedule for step size

    for t = 1:T
        yz = A*z;
        grad  = (1/m)* A'*( ( abs(yz).^2-y ) .* yz ); % Wirtinger gradient
        z = z - mu(t)/normest^2 * grad;             % Gradient update 
    end

    recoveredSoln = z;
    etime = toc;


end

return

