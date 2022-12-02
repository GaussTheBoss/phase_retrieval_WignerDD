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

%% constructMask.m
%
% This function generates the block-circulant (correlation-based) masks 
% used for phase retrieval in 
%
%   Fast Phase Retrieval for High-Dimensions
%   Mark Iwen, Aditya Viswanathan, Yang Wang
%   Jan 2015
%   arXiv:1501.02377
%
% The masks are based on the measurement model
%
%   b^i = |corr(y,w)|^2 + n_i,  i = 1,...,P 
%
% where
%
%   corr(x,y) is the discrete (cross) correlation between x and y
%   w^i \in C^d is a mask having compact support
%   b^i \in R^d are the magnitude measurements
%   n_i \in R^d denotes measurement noise
%   
% The masks have compact support; in particular each mask has only delta+1
% non-zero entries. For a description of mask construction, see the 
% arXiv preprint.
%
% Usage:
%   [maskEntries, blockDiag, measurementMat] = contructMask(pars)
%
% Inputs:
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
%   maskEntries     -   Non-zero mask entries
%   blockDiag       -   Block-diagonalization of this construction
%   measurementMat  -   Full measurement matrix (for use with PhaseLift) 
%

function [maskEntries, blockDiag, measurementMat] = constructMask_EigPR(pars)

%% Initialization

maskType    = pars.maskType;           % Type of mask to generate
if( isfield(pars, 'delta') )
    delta       = pars.delta;              % Block parameter
end
if( isfield(pars, 'over') )
    over        = pars.over;               % Oversampling factor
else
    over = 1;
end
d           = pars.d;                  % Problem size

% No. of masks (this is >= to 2delta+1)
nMasks = round( over*(2*delta-1) );


%% Generate masks

switch lower(maskType)
    case 'random'
        % Random Mask
        % Mask entries are uniform random
        maskEntries = ( randn(delta, nMasks) - 0.5 ) + ...
                                    1i*( randn(delta, nMasks) - 0.5 );
                                
        % The "full" measurement matrix
        if( strcmp(pars.alg, 'PhaseLift') || pars.fullMat)
            % For PhaseLift, store the full matrix
            % Note that PhaseLift is generally used with small problem 
            % sizes, for which using the full matrix is easier and faster 
            % than using FFTs.
            measurementMat = zeros(nMasks*d, d);
            maskEntries = maskEntries.';
            entries = zeros(1,d);
            for ix = 1:nMasks
                entries(1:delta) = maskEntries(ix,:);
                measurementMat( (ix-1)*d+1:ix*d, : ) = ...
                                            gallery( 'circul', entries );
            end
            measurementMat = sparse(measurementMat);
            maskEntries = maskEntries.';
        else
            measurementMat = 0;
        end
                                
                                
    case 'rayan'
        % Rayan/Bryan construction with linear condition number
        maskEntries = zeros(delta, nMasks);
        
        % Standard basis
        ej = zeros(delta,1);
        
        % First mask
        e1 = ej; e1(1) = 1;
        maskEntries(:,1) = e1;
        
        for ix = 1:delta-1            
            ej(ix+1) = 1;             % j^th basis vector
            
            % 
            maskEntries(:,2*ix)   = e1 +    ej;
            maskEntries(:,2*ix+1) = e1 - 1i*ej;
            
            ej(ix+1) = 0;             % reset j^th basis vector    
        end
                
        % The "full" measurement matrix
        if( strcmp(pars.alg, 'PhaseLift') || pars.fullMat )
            % For PhaseLift, store the full matrix
            % Note that PhaseLift is generally used with small problem 
            % sizes, for which using the full matrix is easier and faster 
            % than using FFTs.
            measurementMat = zeros(nMasks*d, d);
            maskEntries = maskEntries.';
            entries = zeros(1,d);
            for ix = 1:nMasks
                entries(1:delta) = maskEntries(ix,:);
                measurementMat( (ix-1)*d+1:ix*d, : ) = ...
                                            gallery( 'circul', entries );
            end
            measurementMat = sparse(measurementMat);
            maskEntries = maskEntries.';
        else
            measurementMat = 0;
        end

        
    case 'fourier'
        % Mark's complex exponential (Fourier-type) masks
        % See description in arXiv preprint
        maskEntries = exp( 2i*pi*(0:delta-1).'*(0:nMasks-1)/nMasks...
            ) .* exp( -(0:delta-1).'*ones(1,nMasks)/delta )/(nMasks^0.25);

        % The "full" measurement matrix
        if( strcmp(pars.alg, 'PhaseLift') || pars.fullMat )
            % For PhaseLift, store the full matrix
            % Note that PhaseLift is generally used with small problem 
            % sizes, for which using the full matrix is easier and faster 
            % than using FFTs.
            measurementMat = zeros(nMasks*d, d);
            maskEntries = maskEntries.';
            entries = zeros(1,d);
            for ix = 1:nMasks
                entries(1:delta) = maskEntries(ix,:);
                measurementMat( (ix-1)*d+1:ix*d, : ) = ...
                                            gallery( 'circul', entries );
            end
            measurementMat = sparse(measurementMat);
            maskEntries = maskEntries.';
        else
            measurementMat = 0;
        end

    case 'global-randgauss'
        % Global random Gaussian measurements
        % No. of measurements is nMasks x d
        D = nMasks*d;
        maskEntries = randn(D, d) + 1i*randn(D, d);
        measurementMat = maskEntries;
        
end
        
% Return the transpose
maskEntries = maskEntries.';
% Structure of the mask entries is as follows:
% maskEntries(i,j) = mask_i, j^th entry
% Note: There are over(2delta-1) masks, each having delta non-zero entries


%% Construct block-circulant matrix for the phase differences
% We will now construct the blocks of the block-circulant matrix M'
%
% Note: Block_j has the form
% [(m_1)j,j (m_1)j,j+1 ... (m_1)j,delta 0 ... 0 (m_1)j+1,1 (m_1)j+1,2 ... (m_1)j+1,j
%  (m_2)j,j (m_2)j,j+1 ... (m_2)j,delta 0 ... 0 (m_2)j+1,1 (m_2)j+1,2 ... (m_2)j+1,j
%                            ....
%  (m_i)j,j (m_i)j,j+1 ... (m_i)j,delta 0 ... 0 (m_i)j+1,1 (m_i)j+1,2 ... (m_i)j+1,j]
% where i denotes the number of masks and (m_i)j,k = (m_i)j (m_i)k^*.
%
% Note: No. of blocks is delta
%       Dimension of each block is over(2delta-1) x (2delta-1)
nMasks = size(maskEntries, 1);
circulantMat_blocks = zeros(nMasks, 2*delta-1, delta);
% Construct each block
% First, the leading entries of a block ...
for ia = 1:delta
    for ib = ia:delta            
        circulantMat_blocks(:,ib-ia+1,ia) = ...
                            maskEntries(:,ia).*conj( maskEntries(:,ib) );
    end
end
% Next, the trailing entries of a block ...
for ia = 1:delta-1
    for ib = 1:ia
        circulantMat_blocks(:,2*delta-1-ia+ib,ia) = ...
                        maskEntries(:,ia+1).*conj( maskEntries(:,ib) );
    end
end


%% Compute block-diagonalization of this matrix

% Since M' is block circulant, it can be block diagonalized; i.e., 
%       U^* M' U = J, where J = diag(J_k) is block diagonal with blocks
%   J_k = \sum_{l=1}^delta M'_l e^{2 pi i k l/d}.
%
% Compute the blocks J_k
blockDiag = zeros(nMasks, 2*delta-1, d);
for ia = 1:d
    blockDiag(:,:,ia) = sum( circulantMat_blocks .* ...
    reshape( ...
        repelem( exp(2i*pi*(ia-1)*(0:delta-1)/d), nMasks, 2*delta-1 ), ...
                        nMasks, 2*delta-1, delta ), 3);
end


end

