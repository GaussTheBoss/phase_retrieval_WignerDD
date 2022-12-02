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

%% generateMeasuremens.m
%
% This function generates (squared) correlation measurements for simulating
% phase retrieval. Measuremens take the form,
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
% For a description of measurement construction, see arXiv preprint
%   Fast Phase Retrieval for High-Dimensions
%   Mark Iwen, Aditya Viswanathan, Yang Wang
%   Jan 2015
%   arXiv:1501.02377
%
% Usage:
%   [measurements, pars] = generateMeasurements(testFnc, maskEnries, pars)
%
% Inputs:
%   testFnc     -   Test function whose measurements are to be generated
%                   (complex vector, double)
%   maskEntries -   Non-zero entries of masks used
%                   (complex matrix, double)
%                   Note: Structure of the mask entries is as follows:
%                   maskEntries(i,j) = mask_i, j^th entry
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
%   measurements    -   generated measurements
%   pars            -   parameter structure with signal power and noise
%                       power included
%

function [measurements, pars] = generateMeasurements_EigPR( ...
                                testFnc, maskEntries, pars )

%% Initialization

d = pars.d;                 % Problem dimension
addNoise = pars.addNoise;   % Noisy simulation?
snr = pars.snr;             % If noisy simulation, added snr (dB)

% No. of masks used
nMasks = size(maskEntries, 1); 

% What type of (real) noise do we add?
% TODO: Move this to parameter structure
% noiseType = 'uniform';              % Uniform random noise
noiseType = 'gaussian';             % Gaussian noise


%% Use flattening?
% Flatten signal? If yes, generate flattening matrices
if(pars.useFlattening)
    % The flattening transform can be written as
    %   y_flat = P W B y, 
    %       where
    %           P is a random permutation matrix
    %           W is a unitary matrix (for example, the DFT matrix)
    %           B is a diagonal matrix with random +/-1 entries

    % First, generate random diagonal matrix B
    binVec = binornd(1, .5, d, 1); binVec(binVec==0) = -1;
    B = spdiags( binVec, 0, d, d );

    % Now, generate permutation matrix P
    P = speye(d);
    idx = randperm(d);
    P = P(idx, :);
    
    % Store fast J-L flattening matrices in the parameter structure for
    % later use 
    pars.B = B;
    pars.P = P;
 
    % Flatten signal
    testFnc = P*fft(B*testFnc);
end

%% Measurements

% For PhaseLift, we have the full measurement matrix...
if( strcmp(pars.alg, 'PhaseLift') || ...
                        strcmp(pars.maskType, 'global-randGauss') )
    
    measurements = abs(maskEntries*testFnc).^2;

else

% Store measurements in a d x P matrix, each column containing the
% measurements corresponding to a distinct mask
measurements = zeros(d,nMasks);

% Since we use correlation measurements of the form 
%   |corr(y,w)|^2 = b^i, i = 1,...,P, 
% these can be efficiently computed using FFTs. This is especially useful
% when the problem dimension is large.
%
% Principle:
%   We express the correlation as a matrix multiplication
%   Since this correlation matrix has a circulant structure, it can be 
%   generated as 
%       corrMat = gallery('circul', maskEntries)
%   Note: Matlab's gallery command generates a circulant matrix given the
%   entries of the mask
%
%   Now, 
%   measurements = abs( corrMat * trueFunction ).^2;
%   The multiplication (corrMat * trueFunction) can be computed efficiently
%   using FFTs
%   We repeat for measurements corresponding to each mask


% First, pre-store FFT of the true function
fftFun = fft(testFnc);

% For each mask ...
for idx = 1:nMasks
    % Evaluate correlations efficienly
    fastMeasurement = conj(fftFun) .* fft( maskEntries(idx,:).', d );
    fastMeasurement = abs( ifft(fastMeasurement) ).^2;      % Squared 
                                                            % measurements
    measurements(1, idx) = fastMeasurement(1);
    measurements(2:end, idx) = fastMeasurement(end:-1:2);        
end

end

%% Add Noise
if( addNoise )        
    % Calculate signal power
    % This is just the norm^2 of the measurements divided by dimension
    signal_power = norm( measurements(:) )^2 / numel(measurements);
    % Store this in the parameters
    pars.signal_power = signal_power;
    
    % Calculate noise power to be added
    noise_power = signal_power / ( 10^(snr/10) );
    % Store this in the parameters
    pars.noise_power = noise_power;
    
    % Generate the noise vector
    switch lower(noiseType)
        case 'uniform'
            % Uniform random noise            
            % Noise model: Additive uniform random noise
            % b_i^ = b_i + n_i, where
            %   b_i is ith entry of measurement vector
            %   n_i is the added noise
            %   b_i^ is the ith entry of noise corrupted measurement
            % n_i ~ U[0,a] is iid uniform random noise
            % Note: mean(n) = a/2, var(n) = (a^2)/12
            % We choose a as follows:
            %   mean^2 + var^2 = noise_power
            %   (1/4 + 1/12)*a^2 = a^2/3 = noise_power

            % Choose uniform distribution parameter a
            alpha = sqrt( 3*noise_power );
            noise = alpha*rand( size(measurements) );
            
        case 'gaussian'
            % Add (real) Gaussian noise of desired variance
            noise = sqrt(noise_power)*randn( size(measurements) );
            
    end    
    
    % Add noise to measurements
    measurements = measurements + noise;    
end


end

