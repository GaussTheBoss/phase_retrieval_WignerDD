function xrec = angsync( X, d, kappa, mags )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% Initial guess
% Rather than start with a random guess, lets start with one of the
% standard basis elements, e_k, k = max_i |x_i|
v_k = zeros(d, 1);
[~, maxloc] = max(mags);
v_k(maxloc) = 1;
    
% Power method - compute top e-vec
% implementation of the shifted inverse power method
maxit_stage1 = 50;              % max iterations - first stage
maxit_stage2 = 50;
tol = 1e-13;                    % tolerance
    
% we know top eivenvalue is 2kappa-1
mu = (2*kappa-1);
X_kappa = X-mu*speye(d);

[L, U, P] = lu(X_kappa);           % precompute LU decomposition
% Shifted power method
for ia = 1:maxit_stage1
    v_prev = v_k;
    w = U\(L\(P*v_k));              % solve (A - mu*I)w = v_k
    v_k = w/norm(w);
    v_k = v_k./abs(v_k);
    % if tolerance condition satisfied, exit loop...
    if( norm(v_prev-v_k) <= tol ) 
        break;  
    end
end
    
% Rayleigh quotient refinement/iteration
for ib = 1:maxit_stage2
    v_prev = v_k;
    % Rayleigh quotient
    mu = v_k'*(X*v_k)/(v_k'*v_k);
    X_kappa = X-mu*speye(d);
    w = X_kappa\v_k;
    v_k = w/norm(w);
    v_k = v_k./abs(v_k);
    % if tolerance condition satisfied, exit loop...
    if( norm(v_prev-v_k) <= tol )
        break;  
    end        
end


nz_idx = find(v_k);       % of only the non-zero entries
v_k(nz_idx) = v_k(nz_idx) ./ abs(v_k(nz_idx));
xrec = v_k;


end

