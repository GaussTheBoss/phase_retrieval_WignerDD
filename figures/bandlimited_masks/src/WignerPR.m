function [xrec, etime] = WignerPR(Y, pars, mask)

%% Initialization
% Read in problem parameters
d = pars.d;
L = pars.L; 
kappa = pars.kappa;
pproc = pars.pproc;

mask_precomp_ld = mask.ld;
mask_precomp_ud = mask.ud;
m = mask.m;

                              
%% Solve for diagonal of x^x^*
tic;
% First, compute the left-hand  side of double-aliasing formulation
% TODO: add Sami's paper reference
LHS = fft( fft(Y,L,2),d ).';

% Next, solve for diagonals
X_diags = zeros(d,2*kappa-1);
for alpha = 0:kappa-1
%     X_diags(:,alpha+1) = ((d^2)/L) * LHS(alpha+1,:).' ./ ...
%                                             mask_precomp_ld(:,alpha+1);  
    X_diags(:,alpha+1) =  LHS(alpha+1,:).' ./ ...
                                            mask_precomp_ld(:,alpha+1);  

end

lb_idx = L-(kappa-1);
for alpha = L-(kappa-1):L-1    
    % TODO: <check indexing for LHS>
%     X_diags(:,kappa+alpha-lb_idx+1) = ((d^2)/L) * LHS(alpha+1,:).' ...
%                                     ./ mask_precomp_ud(:,alpha-lb_idx+1);
    X_diags(:,kappa+alpha-lb_idx+1) = LHS(alpha+1,:).' ...
                                    ./ mask_precomp_ud(:,alpha-lb_idx+1);
                                
end

X_diags = conj(ifft(X_diags, d));

% Put into diagonal form
diagLocs = [0:kappa-1 -(kappa-1):-1 ...
            -(d-1):-(d-1)+(kappa-2) (d-1)-(kappa-2):(d-1)];
X = spdiags([X_diags X_diags(:,2:end)], diagLocs, d, d);


% Hermitian symmetrize
% TOD): <fix the conj(.)>
X = X/2 + X'/2;

% Magnitudes
% if (pars.magEst)
switch lower(pars.magEst)
    case 'eig'
        mags = (d/sqrt(L))*EstimateMagnitude( X, pars );
    case 'diag'
        mags = (d/sqrt(L))*sqrt( diag(X) );
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
        Z = (d/sqrt(L))*(X - diag(diag(X)));
        W = abs(Z);         % weights on graph
        gL = diag(W*ones(d,1))-W.*sign(Z);   % graph Laplacian

        % this is Matlab's eig solver
        [xrec, ~, ~] = eigs(gL, 1, 'SM');
        xrec = sign(xrec);
        
end

xrec = ifft( full(mags.*xrec), d );


%% Post-Processing (HIO+ER)

etime = toc;

if( pproc)
    [xrec, pproctime] = postprocess_HIOER( xrec, Y, m, d, L );
    etime = etime + pproctime;                                    
end
    
return