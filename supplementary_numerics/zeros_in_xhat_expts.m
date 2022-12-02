% Script to check robustness of angular synchronization to presence of
% zeros in x^hat

clear; close all; clc
addpath ./src/

% For repeatability, set random seed to ...
rng(12345);


%% Parameters
d =  2^08;                      % signal dimension
rho = round(1.25*log2(d));    % support of m^hat
kappa = round(rho/1);         % Note: kappa should obey L=delta-1+kappa
L = rho+kappa-1;              % subsampling in space - no. of shifts

% Fix d to satisfy divisibility requirements
d = L*round(d/L);

                
% Percent of zeros in xhat
percent_zeros = (0:5:95).';

% No. of trials
ntrials = 100;

% Store error here
error = zeros( length(percent_zeros), 2 );


for izeros = 1:length(percent_zeros)
    
    % Experiments with x% of zeros in xhat
    pzeros = percent_zeros(izeros);
    fprintf( '\n\n Now running tests with %d percent zeros...', pzeros );

    
for itrials = 1:ntrials    

    %% Generate test x^hat with a certain number of zeros
    xhat = randn(d,1)+1i*randn(d,1);                % random complex signal x
    
    % entries randomly zeroed out
    xhat( randperm(d,round(d*pzeros/100)) ) = 0;
        
    % Generate rank-1 matrix/outer product
    X_true = xhat*xhat';

    % ... and zero out all but the required diagonals
    diagLocs = [0:kappa-1 -(kappa-1):-1 ...
                -(d-1):-(d-1)+(kappa-2) (d-1)-(kappa-2):(d-1)];

    X = zeros(d);        
    for ix = 1:length(diagLocs)
        entries = diag(X_true, diagLocs(ix));
        X = X + diag( entries, diagLocs(ix) );
    end

    % Hermitian symmetrize
    X = X/2 + X'/2;
    
    % We will use diagonal magnitude estimation and use eigenvector based 
    % angular synchronization
    magEst = 'diag';
    angSync = 'topeig';

    % Magnitudes
    switch lower(magEst)
        case 'eig'
            pars.d = d; pars.kappa = kappa;
            mags = EstimateMagnitude( X, pars );
        case 'diag'
            mags = sqrt( diag(X) );
    end



    %% Angular synchronization
    switch lower(angSync)
        case 'topeig'
            % entrywise normalize
            nz_idx = find(X);
            X(nz_idx) = X(nz_idx)./abs(X(nz_idx));

            % % this uses the shifted inverse power method
            % xrec = angsync( X, d, kappa, mags );

            % this is Matlab's eig solver
            [xrec, ~, ~] = eigs(X, 1, 'LM');
            xrec = sign(xrec);

        case 'graphleig'
            % find eigenvector corr. to smallest eigenvalue of graph Laplacian
            Z = (d/sqrt(L))*(X - diag(diag(X)));
            W = abs(Z);         % weights on graph
            gL = diag(W*ones(d,1))-W.*sign(Z);   % graph Laplacian

            % this is Matlab's eig solver
            [xrec, ~, ~] = eigs(gL, 1, 'SM');
            xrec = sign(xrec);

    end


    % Recovered signal after ang. sync.
    xhatrec = mags.*xrec;

    % Correct for global phase factor
    phaseOffset = angle( (xhatrec'*xhat) / (xhat'*xhat) );
    xhatrec = xhatrec * exp(1i*phaseOffset);

    % Reconstruction error
    errordB = 10*log10( norm(xhat-xhatrec)^2 / norm(xhat)^2 );
    relerror = ( norm(xhat-xhatrec) / norm(xhat) );
    error(izeros, 1) = error(izeros, 1) + errordB/ntrials;
    error(izeros, 2) = error(izeros, 2) + relerror/ntrials;
    
end

    fprintf( '\n  (Avg.) error in Ang. Sync. is %3.2fdB', error(izeros, 1) );
    fprintf( '\n  (Avg.) error in Ang. Sync. is %3.2e\n', error(izeros, 2) );

end



% plot
figure; 
semilogy(percent_zeros, error(:,2), '-d', 'linewidth', 2);

axis([-5 105 1e-15 1e1]); grid
% xlabel( 'Percent of Zero Entries in Signal', 'interpreter', 'latex', 'fontsize', 14 )
xlabel( 'No. of Consecutive Zero Entries in $\widehat{\mathbf{x}}$', 'interpreter', 'latex', 'fontsize', 14 )
ylabel( 'Relative 2-norm Error', 'interpreter', ...
                                                'latex', 'fontsize', 14 )
title( 'Error in Angular Synchronization, $d=247, \rho=\kappa=10$', ...
                                 'interpreter', 'latex', 'fontsize', 14 )

% Save data
% save( 'zeros_in_xhat_expt.mat', 'percent_zeros', 'error' );