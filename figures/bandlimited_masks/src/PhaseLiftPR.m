function [xrec, etime] = PhaseLiftPR(Y, pars, measurementMat, varargin)

%% Initialization
% Read in problem parameters
d = pars.d;
L = pars.L;

% Process optional arguments
if(nargin == 4)
    tfocsLinOp = varargin{1};
end

                              
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


%--------------------------------------------------------------------------
   % We'll implement this as a trace regularized least-squares problem in
    % TFOCS
    
    lambda = 1e-6;
    % Set TFOCS parameters
    % TODO: Move this to main parameters structure
    opts.maxIts     = 1e3;
    opts.tol        = 1e-10;
    opts.restart    = 200;
    % opts.largescale = true;
    opts.printEvery = 0;        % Suppress TFOCS status messages/output
    
    % Initial guess
    initGuess = zeros(d);
    
%     % TFOCS requires special initialization of the measurement model
    if( nargin~=4 )
        tfocsLinOp = initializeLinopPR( {[d, d], [L*d, 1]}, measurementMat ); 
    end
    X = solver_TraceLS( tfocsLinOp, Y, lambda, initGuess, opts );

%--------------------------------------------------------------------------


% The above SDP recovers the matrix X = xx*; we need to extract x
% Since X is low-rank (ideally rank-1), choose solution to be (scaled) 
% leading eigenvector                                        
[recoveredSoln, eVal] = eig(X);
xrec = sqrt(eVal(end))*recoveredSoln(:,end);
    
etime = toc;                                    

return