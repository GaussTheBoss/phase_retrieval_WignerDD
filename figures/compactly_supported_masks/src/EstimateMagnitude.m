% Copyright (c) 2016- Michigan State University and the CHARMS Research
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

%% EstimateMagnitude.m
%
% Magnitude estimation from scaled phase difference measurements.
%
% Let x \in C^d and define X \in C^{d x d} to be the (possibly noise 
% corrupted) banded matrix
%
%                   {   (xx*)j,k        if |j- k mod d| < \delta 
%       (X)j,k =    {
%                   {       0           else
%
% This function estimates |x_i|, i=1,...,d. The following methods are
% available:
%       (1)     Diagonal estimation
%                   Sets |x_i| = sqrt( |(X)j,j| )
%       (2)     Eigen estimation
%                   Computes the leading eigenvector of the \delta x \delta
%                   blocks
%
%          [    (X)j,j       (X)j,j+1       .....    (X)j,j+\delta     ]
%          [   (X)j+1,j     (X)j+1,j+1      .....   (X)j+1,j+\delta    ]
%          [       .             .            .             .          ]
%          [       .             .            .             .          ]
%          [       .             .            .             .          ]
%          [ (X)j+\delta,j (X)j+\delta,j+1  ..... (X)j+\delta,j+\delta ]
%
%                   for j = 1,....,d and averages the result. This gives a
%                   more noise robust estimate of the magnitudes. Also note
%                   that the block size could be between 1 and delta and j
%                   could be advanced in increments of between 1 and delta.
%
%
% Usage:
%   EstMag = EstimateMagnitude(bandedMat, pars)
%
% Inputs:
%   bandedMat   -   Banded matrix of scaled Phase differences. See above
%                   for structure. Each row has 2\delta-1 non-zero entries 
%                   of the form (X)j,k = x_jx_k^*.
%                   Matrix, double, complex, dxd
%
%   pars        -   Matlab structure containing problem parameters.
%                   See <root>/src/PrintProblemPars.m for elements.
%
% Output(s):
%   EstMag      -   Estimated magnitude vector.
%                   Vector, double, real, dx1
%

function EstMag = EstimateMagnitude( bandedMat, pars)

%% Initialization

d           = pars.d;                       % Problem dimension
kappa       = pars.kappa;                   % Block/Band parameter
MagEstType = 'eig';                         % Default mag. estimation

%% Magnitude Estimation
switch lower(MagEstType)
    
    %% Diagonal magnitude estimation
    case 'diag'
        EstMag = sqrt( abs(diag(bandedMat)) );
    
    %% Eigenvector-based magnitude estimation
    case 'eig'
        % Parameters
        blkLength   = kappa;    % Process matrix in blocks of this size
        shift       = 1;        % Overlapping block shift value
        
        % Matlab eigenvalue solver (eigs) options
%         opts.isreal = false; opts.issym = true;
    %     opts.tol = 1e-6;
    %     opts.v0 = zeros(blkLength, 1);
        % Circular indexing
        wrap = @(x) mod(x-1, d) + 1;

        % Store magnitude estimate here
        EstMag = zeros(d, 1);
        
%         % Use diagonal estimate to initilize eigenvector method
%         DiagEst = sqrt( abs(diag(bandedMat)) );
        
        for ia = 1:shift:d
            % Extract block
            blkMat = bandedMat( wrap(ia:ia+blkLength-1), ...
                                                wrap(ia:ia+blkLength-1) );
            % Compute top eigenvector
%             opts.v0 = DiagEst( wrap(ia:ia+blkLength-1) );
            [evec, eval] = eigs(blkMat, 1, 'LM');
        
            % Update magnitude estimate
            % This averages the estimates from multiple overlapping blocks
            EstMag( wrap(ia:ia+blkLength-1) ) = ...
                EstMag( wrap(ia:ia+blkLength-1) ) + ...
                ( sqrt(abs(eval)) * abs(evec) )/ceil(blkLength/shift);

        end
        
end


return