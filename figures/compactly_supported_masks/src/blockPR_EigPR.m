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

function [recoveredSoln, execTime] = blockPR_EigPR( measurements, ...
                                            maskEntries, blockDiag, pars)

%% blockPR.m
%
% For details, see arXiv preprint
%   Fast Phase Retrieval for High-Dimensions
%   Mark Iwen, Aditya Viswanathan, Yang Wang
%   Jan 2015
%   arXiv:1501.02377
%
% Solve the phase retrieval problem using (squared magnitude) local 
% correlation-based phaseless measurements; i.e., 
%
%   find  x    from     |M x|^2 = b, 
%
% where
%   x \in C^d is the unknown signal
%   b \in R^d are the magnitude measurements
%   M \in C^{Dxd} is block-circulant, with 
%                [     M_1     ]
%                [     M_2     ]
%         M =    [      .      ]
%                [      .      ]
%                [      .      ]
%                [ M_{2delta-1}]
%    where 
%       M_i \in C^{dxd} are circulant matrices describing correlations 
%       with local (compactly-supported) masks. 
%
%       For example, let d=4 and let each measurement mask have only 
%       delta=2 non-zero entries. Then
%           m_i \in C^4 = [ (m_i)_1   (m_i)_2    0    0]^T
%       and
%                    [ (m_1)_1^*  (m_1)_2^*      0          0     ]
%             M_i =  [     0      (m_1)_1^*  (m_1)_2^*      0     ]
%                    [     0           0     (m_1)_1^*  (m_1)_2^* ]
%                    [ (m_1)_2^*       0          0     (m_1)_1^* ]
%       Total no. of measurements is D = (2delta-1)d = 12.
%
%
% Usage:
%   [recoveredSoln, execTime, precompTime] = blockPR( trueFunction, ...
%                                       measurements, pars)
%
% Inputs:
%   trueFunction    -  True function (complex vector, double)
%   measurements    -  Acquired measurements (real array, double)
%                       array is of dimension d x nMasks, where 
%                       nMasks is the no. of masks used (this depends on
%                       the block parameter delta and oversampling factor).
%   blockDiag       -  Pre-computed block-diagonalization of the
%                      measurement construction
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
%   recoveredSoln   -  Recovered solution (complex vector, double)
%   preCompTime     -  Pre-computation time
%   execTime        -  BLockPR execution time
%

%% Initialization

d = pars.d;                                 % Problem size
delta = pars.delta;                         % Band parameter, delta

% No. of masks used
nMasks = round( pars.over*(2*delta-1) );


               
%% BlockPR - Phase Retrieval from Block-Circulant Measurements
%
% The algorithm can be conceptually divided into two steps:
%   (i) Banded lifting: Recover phase differences (on a band) by solving a 
%       linear system. This returns a banded version of the rank-1 matrix
%       xx^*
%   (ii) Angular Synhronization: Recover the vector x from the banded phase
%       difference matrix
%

%% Banded Lifting
%       Banded lifting: Recover phase differences (on a band) by solving a 
%       linear system. This returns a banded version of the rank-1 matrix
%       xx^*
%
% See preprint on arXiv for more details. Below is a summary of the system
% we solve for delta=2, d=4
%       M'z = b, where 
%   
% TODO: Write down M', z and b
%

% Actual solution starts here; so record time
tic;

% The block circulant structure of the matrix means that we can efficiently
% invert the system matrix using FFTs.
% (see preprint for proof that system matrix is well conditioned)

%   (a) Begin by processing RHS
%       We take the FFT of each measurement (corr. to a single mask) and 
%       interleave/stack the measurements together
b_hat = fft(measurements)/sqrt(d);
b_hat = reshape( b_hat.', numel(measurements), 1 );

%   (b) Solve a block diagonal system
%       ... with blocks J_k
x_hat = zeros((2*delta-1)*d, 1);
for ia = 1:d
    x_hat( (ia-1)*(2*delta-1)+1:ia*(2*delta-1) ) = blockDiag(:,:,ia) \ ...
                        b_hat( (ia-1)*nMasks+1:ia*nMasks );
end

%   (c) Multiply (using FFTs) with the unitary matrix, U
% First, break up the long vector into vectors corresponding to each mask
x_hat = reshape( x_hat, (2*delta-1), d );
% Compute the inverse FFT
phaseDiffs = sqrt(d)*ifft(x_hat, [], 2)';

% Reshape/unwrap phase differences
for ia = 1:delta-1
    phaseDiffs(:,ia+1) = circshift(phaseDiffs(:,ia+1), ia);
    phaseDiffs(:,end-ia+1) = circshift(phaseDiffs(:,end-ia+1), -(ia-1));
end


%% Angular Synchronization
%       Angular Synhronization: Recover the vector x from the banded phase
%       difference matrix
recoveredSoln = angularSync_EigPR( phaseDiffs, pars, measurements, maskEntries );


%% Post-processing

if( pars.altProj )
%     % Alternating Projection post-processing
%     recoveredSoln = altProj_BlockPR( ...
%                     recoveredSoln, ...      % Initial guess
%                     measurements, ...       % Measurements
%                     maskEntries, ...        % Measurement masks
%                     pars ...                % Problem parameters
%                     );

    % HIO+ER post-processing
    recoveredSoln = HIO_BlockPR_EigPR( ...
                    recoveredSoln, ...      % Initial guess
                    measurements, ...       % Measurements
                    maskEntries, ...        % Measurement masks
                    pars ...                % Problem parameters
                    );

end

% % Revert flattening transform (this is typically used with sparse signals)
% if( pars.useFlattening )
%     % The inverse flattening transform can be written as
%     %   y = B' W P' y_flat, 
%     %       where
%     %           P is the random permutation matrix
%     %           W is the unitary matrix (for example, the DFT matrix)
%     %           B is the diagonal matrix with random +/-1 entries
%     recoveredSoln = pars.B'*ifft( full(pars.P'*recoveredSoln) );
% end

% Record execution time
execTime = toc;

return