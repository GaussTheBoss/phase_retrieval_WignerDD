function [xrec, etime] = WignerPR(Y, pars, mask, varargin)

%% Initialization

if(nargin > 3)
    maskpars = varargin{1};
    maskEntries = varargin{2};
end

% Read in problem parameters
d = pars.d;
K = pars.K; 
kappa = pars.kappa;
pproc = pars.pproc;

mask_precomp_ld = mask.ld;
mask_precomp_ud = mask.ud;
m = mask.m;

tic; 


% First, compute the left-hand  side of double-aliasing formulation
% TODO: add Sami's paper reference
LHS = fft( fft(Y,d,2), K ).'/K;

% Next, solve for diagonals
X_diags = zeros(d,2*kappa-1);
for omega = 0:kappa-1
    X_diags(:,omega+1) =  LHS(:,omega+1) ./ mask_precomp_ud(:,omega+1);  

end

lb_idx = K-(kappa-1);
for omega = K-(kappa-1):K-1    
    % TODO: <check indexing for LHS>
    X_diags(:,kappa+omega-lb_idx+1) = LHS(:,omega+1) ...    
                                    ./ mask_precomp_ld(:,omega-lb_idx+1);
                                
end


X_diags = (ifft(X_diags, d));


% Put into diagonal form
X = spdiags(X_diags(:,1), 0, d, d);
% first the upper diagonals
for ix = 1:kappa-1
    X = spdiags(circshift(X_diags(:, ix+1),ix), ix, X);
    % this is the circulant bits
    X = spdiags(circshift(X_diags(:, ix+1), -(d-ix)), -(d-ix), X);
end

% ... and now the lower diagonals
for ix = 1:kappa-1
    X = spdiags(circshift(X_diags(:, end-ix+1),-ix), -ix, X);
    % this is the circulant bits
    X = spdiags(circshift(X_diags(:, end-ix+1), (d-ix)), (d-ix), X);
end


% Hermitian symmetrize
X = X/2 + X'/2;

% Magnitudes
% if (pars.magEst)
switch lower(pars.magEst)
    case 'eig'
        mags = EstimateMagnitude( X, pars );
    case 'diag'
        mags = sqrt( diag(X) );
end



%% Angular synchronization
switch lower(pars.angSync)
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
        Z = X - diag(diag(X));
        W = abs(Z);                         % weights on graph
        gL = diag(W*ones(d,1))-W.*sign(Z);  % graph Laplacian

        % this is Matlab's eig solver
        [xrec, ~, ~] = eigs(gL, 1, 'SM');
        xrec = sign(xrec);
        
end

xrec = xrec.*mags;


%% Post-Processing (HIO+ER)
    
etime = toc;

if( pproc)
    
%     xrec = postprocess_HIOER( xrec, Y, m, d, L, [120, 25, 5] );
%     xrec = postprocess_HIOER( xrec, Y, m, d, K );
%     xrec = HIO_BlockPR( xrec, Y.', maskEntries, maskpars );    
    
    [xrec, pproctime] = postprocess_HIOER( xrec, Y, m, d, K );
    etime = etime + pproctime;                                    
end


return