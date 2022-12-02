function [xrec, etime] = PhaseLiftPR(Y, pars, tfocsLinOp)

%% Initialization
% Read in problem parameters
d = pars.d;

                              
%% Here is the PhaseLift solution
% We implement PhaseLift as a trace regularized least-squares problem in
% CVX

% % Reg. parameter
% lambda = 1e-4;
% 
% cvx_begin sdp
%     % If you do not have access to the MOSEK solver (MOSEK ApS or CVX 
%     % Professional license), comment this line
%     cvx_solver mosek
%     cvx_quiet true;
%     variable X(d,d) hermitian
%         
%     minimize 0.5*norm( diag(measurementMat * X * measurementMat') - Y(:) )...
%                         + lambda*trace( X )
%     subject to
%         X >= 0;
% cvx_end

tic;


%--------------------------------------------------------------------------
   % We'll implement this as a trace regularized least-squares problem in
    % TFOCS
    lambda = 1e-6;
%    switch pars.snr
%        case 10
%            lambda = 1e-6;
%            opts.maxIts     = 1.5e3;
%        case 20
%            lambda = 1e-6;
%            opts.maxIts     = 1e3;
%        case 30
%            lambda = 1e-6;
%            opts.maxIts     = 2e3;
%        case 40
%            lambda = 1e-6;
%            opts.maxIts     = 3e3;
%        case 50
%            lambda = 1e-4;
%            opts.maxIts     = 4e3;
%        case 60
%            lambda = 1e-6;
%            opts.maxIts     = 5e3;
%    end

   switch pars.snr
       case 10
           lambda = 1e-6;
           opts.maxIts     = 1e3;
       case 20
           lambda = 1e-6;
           opts.maxIts     = 1e3;
       case 30
           lambda = 1e-6;
           opts.maxIts     = 1.2e3;
       case 40
           lambda = 1e-6;
           opts.maxIts     = 1.9e3;
       case 50
           lambda = 1e-4;
           opts.maxIts     = 2e3;
       case 60
           lambda = 1e-6;
           opts.maxIts     = 5e3;
   end


    % Set TFOCS parameters
    % TODO: Move this to main parameters structure
%     opts.maxIts     = 1e3;
    opts.tol        = 1e-12;
    opts.restart    = 200;
    opts.largescale = true;
    opts.printEvery = 0;        % Suppress TFOCS status messages/output
    
    % Initial guess
    initGuess = zeros(d);
    
%     % TFOCS requires special initialization of the measurement model
%     tfocsLinOp = initializeLinopPR( {[d, d], [K*d, 1]}, ...
%                                                     measurementMat ); 
    X = solver_TraceLS( tfocsLinOp, Y, lambda, initGuess, opts );

%--------------------------------------------------------------------------


% The above SDP recovers the matrix X = xx*; we need to extract x
% Since X is low-rank (ideally rank-1), choose solution to be (scaled) 
% leading eigenvector                                        
[recoveredSoln, eVal] = eig(X);
xrec = sqrt(eVal(end))*recoveredSoln(:,end);
    
etime = toc;                                    

return
