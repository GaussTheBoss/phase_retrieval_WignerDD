% Copyright (c) 2015 Michigan State University and the CHARMS Research
% Group
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

%% angularSync.m
%
% This function solves the angular synchronization problem,
%
% Estimate d unknown angles theta_1, theta_2, ... , \theta_d \in [0, 2π) 
% from d(2\delta - 1) (noisy) measurements of their differences
%   \Delta \theta_{ij} := \theta_i - \theta_j
%
%
% Usage:
%   recoveredSoln = angularSync( phaseDiffs, pars)
%
% Inputs:
%   phaseDiffs -    Phase differences
%                   TODO: Summarize structure of this array
%   pars    -   Matlab structure containing problem parameters.
%               Structure elements are
%                   'd'             - Problem dimension (integer)
%                   'delta'         - Block length/parameter (integer)
%                   'addNoise'      - Noisy simulation? (boolean)
%                   'snr'           - If noisy simulation, added noise SNR 
%                                     (in dB) (double)
%                   'testFncType'   - type of test function (string)
%                   'nComps'        - for particular test functions 
%                                     (integer)
%                   'useFlattening' - use fast JL flattening? (boolean)
%                   'maskType'      - type of masks (string)
%                   'over'          - oversampling factor (double)
%                   'angSyncMethod' - angular synchronization method
%                                     (string)
%                   'threshold'     - noise threshold for angular
%                                     synchronization algorithms (double)
%                   'altProj'       - use alternating projection
%                                     post-processing? (boolean)
%
% Output(s):
%   recoveredSoln   -   recovered solution (complex vector, double)
%

function recoveredSoln = angularSync_EigPR( phaseDiffs, pars, varargin)

%% Initialization

d = pars.d;                     % Problem dimension
delta = pars.delta;             % Block parameter

% Process optional arguments
if(nargin >=3)
    measurements = varargin{1};
end
if(nargin >=4)
    maskEntries = varargin{2};
end

% Reshape given phase diffs. to obtain a banded version of the matrix 
% xx^*
diagLocs = [0:delta-1 -(delta-1):-1 ...
            -(d-1):-(d-1)+(delta-2) (d-1)-(delta-2):(d-1)];
phaseDiffs_banded = spdiags([phaseDiffs phaseDiffs(:,2:end)], ...
                                            diagLocs, d, d);

% Hermitian Symmetrize
phaseDiffs_banded = (phaseDiffs_banded + phaseDiffs_banded')/2;

% Record magnitudes
% magnitudes = sqrt( abs(real(diag(phaseDiffs_banded))) );

if( ~isfield(pars, 'MagEstType') )
    pars.MagEstType = 'eig';
    % pars.MagEstType = 'diag';
end
magnitudes = EstimateMagnitude_EigPR( phaseDiffs_banded, pars);


% Normalize entries to unit magnitude
nz_idx = find(phaseDiffs_banded);
% for thresholding purposes, record magnitudes ...
idx_threshold = find( abs(phaseDiffs_banded(nz_idx)) <= pars.threshold );

phaseDiffs_banded(nz_idx) = phaseDiffs_banded(nz_idx) ./ ...
                                abs( phaseDiffs_banded(nz_idx) );
phaseDiffs_banded( idx_threshold ) = 0;
    

switch lower(pars.angSyncMethod)
    
    %% Eigenvector-based angular synchronization
    case 'eig'

    % Options to Matlab's eigen solver
    opts.isreal = false;            % We have complex entries in our matrix
%     opts.issym = true;              % Hermitian symmetric
    %TODO: Is this the best number of iterations?
%     opts.maxit = 500;
%     opts.maxit = max(20, floor(10.0*d/(2*delta)) );   % Max. iterations
%     if( pars.addNoise )             % In noisy simulations, only run eigen 
        % TOFIX                     % solver to roughly the level of noise
%         opts.tol = 1e-3*sqrt( 1/( 10^(pars.snr/10) ) );
%     end    
        
%     % Initial guess
%     % Rather than start with a random guess, lets start with one of the
%     % standard basis elements, e_k, k = max_i |x_i|
%     e_k = zeros(d, 1);
%     [~, maxloc] = max(magnitudes);
%     e_k(maxloc) = 1;
%     opts.v0 = e_k;
%     
%     % Power method - compute top e-vec
%     % implementation of the shifted inverse power method
%     v_k = opts.v0;                  % starting guess
%     maxit_stage1 = 50;              % max iterations - first stage
%     tol = 1e-12;                    % tolerance
%     if(pars.addNoise)
%         tol = 1e-3*pars.noise_power;
%     end
%     
%     % we know top eivenvalue is 2delta-1
%     mu = (2*delta-1);
%     A_delta = phaseDiffs_banded-mu*speye(d);
%         
%     [L, U, P] = lu(A_delta);           % precompute LU decomposition
%     % Shifted power method
%     for ia = 1:maxit_stage1
%         v_prev = v_k;
%         w = U\(L\(P*v_k));              % solve (A - mu*I)w = v_k
%         v_k = w/norm(w);
%         % if tolerance condition satisfied, exit loop...
%         if( norm(v_prev-v_k) <= tol ) 
%             break;  
%         end
%     end
%     
%     % Rayleigh quotient refinement/iteration
%     maxit_stage2 = 50;
%     for ib = 1:maxit_stage2
%         v_prev = v_k;
%         % Rayleigh quotient
%         mu = v_k'*(phaseDiffs_banded*v_k)/(v_k'*v_k);
%         A_delta = phaseDiffs_banded-mu*speye(d);
%         w = A_delta\v_k;
%         v_k = w/norm(w);
%         % if tolerance condition satisfied, exit loop...
%         if( norm(v_prev-v_k) <= tol )
%             break;  
%         end        
%     end
% %     [ia ib]
%     
%     % final e-vec solution
%     recoveredSoln = v_k;
    
    % Use Matlab's implementation
    [recoveredSoln, eval, flag] = eigs(phaseDiffs_banded, 1, 'LM', opts);    
    if(flag ~= 0 )
        fprintf( '\n **Error!** eigs did not converge! \n' );
    end


    % Retain the phase of the eigen vector solution
    nz_idx = find(recoveredSoln);       % of only the non-zero entries
    recoveredSoln(nz_idx) = recoveredSoln(nz_idx) ./ ...
                                        abs(recoveredSoln(nz_idx));
    recoveredSoln = recoveredSoln .* magnitudes;
%     recoveredSoln = full(recoveredSoln);
   
    
    %% Greedy angular synchronization
    case 'greedy'

    % Initialize solution (we will fill in the phase angle here)
    recoveredSoln = zeros(d,1);

    % Circular indexing function
    wrap = @(x) mod(x-1,d)+1;

    % Summary of the framework
    % We start with the row corresponding to the largest magnitude entry 
    % The phases of the next delta successive components can be fixed by 
    % looking at the phases of the entries along this row of the matrix
    %
    % Among these delta components, we pick the one with the largest magnitude.
    % This will decide the phase of the next set of components. The process
    % repeats until no more phases need to be set.

    % TODO: Right now, this only considers phase differences x_ix_j^*, j>i.
   
    
    % to ensure we do not overwrite already set phase angles
    phaseFlag = ones(d, 1);
    
    % Start with the largest magnitude element
    mag_copy = full(magnitudes);          % Work with a copy
    [~, loc] = max(mag_copy);
    recoveredSoln(loc) = 0;
    phaseFlag(loc) = 0;
    mag_copy(loc) = 0;
    
    while( nnz(phaseFlag) > 0 )

        % This is the delta-block of entries whose phase is to be set
        blockEntries = wrap( [loc-(delta-1):loc-1 loc+1:loc+(delta-1)].' );
        
        % Now use the phase difference estimates from xx^* to set the phase
        % But we don't want to overwrite already set phases
        recoveredSoln(blockEntries) = ...
           ( 1-phaseFlag(blockEntries) ).*recoveredSoln(blockEntries) + ...
           phaseFlag(blockEntries) .* ( recoveredSoln(loc) - ...
                        angle( phaseDiffs_banded(loc,blockEntries) ).' );
        phaseFlag(blockEntries) =  0;
        
        % Update pointers
        [~, maxloc] = max( mag_copy(blockEntries) );
        loc = blockEntries(maxloc);
        mag_copy(loc) = 0;

    end 

%     % Threshold small entries since they can lead to erroneous phase estimates
%     magnitudes( abs(magnitudes)<=pars.threshold ) = 0;

    % Estimate of the recovered vector
    recoveredSoln = magnitudes .* exp(1i*recoveredSoln);    

end

return