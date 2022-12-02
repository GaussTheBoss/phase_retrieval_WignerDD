% Copyright (c) 2014-2015 Michigan State University and the CHARMS Research
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

function projSoln = HIO_BlockPR( initGuess, measurements, maskEntries, pars )


%% Initialization
%if( pars.addNoise )
%    pars.threshold = 1e-3*pars.noise_power;
%end

switch lower(pars.maskType)
    
    case {'fourier', 'random', 'rayan'}
        % These are local masks
        %% Initializations

        % Problem dimension
        d = pars.d;

        % No. of masks used
        nMasks = size(maskEntries, 1);

        % We will work with absolute value (and not squared absolute value) of the
        % measurements
        mag_meas = sqrt(measurements);

%         % Let threshold vary with noise level
%         pars.threshold = pars.noise_power*1e-1;
        
        % Store result of applying forward operator/measurement matrix here
        y = zeros(d, nMasks);
        P_RA = zeros(d, nMasks);

        % Store fft of the masks
        maskFFT = zeros(d, nMasks);
        for idMasks = 1:nMasks
            maskFFT(:,idMasks) = fft( conj(maskEntries(idMasks, :).'), d );
        end

        % Normalization factor for aplpying the (pseudo) inverse transform
        normFac = sum( abs(maskFFT).^2, 2 );

        % Max. iterations
        nIter = pars.altProjIters;  % assuming this is a multiple of 30
        % Default no. of HIO vs ER iterations
        nHIO = 25; nER = 5;        % 20 iter. of HIO + 10 iter. of ER 
        % If this is otherwise specified ...
        if isfield(pars, 'niter_HIO')
            nHIO = pars.niter_HIO;
        end
        if isfield(pars, 'niter_ER')
            nER = pars.niter_ER;
        end
        nPasses = round(nIter/(nHIO+nER));
        
        % Initial guess
        x = full(initGuess);

        % First, apply forward operator
        %% Forward operator
        % For efficiency, compute these using FFTs
        % Recall that our measurements are correlations, which can be
        % efficiently evaluated using FFTs
        % Pre-store FFT of the current guess
        fftFun = ifft(x);

        for idMasks = 1:nMasks
            % Measurements 
            % Structure of measurements - |corr(y,w)| = b^i, i = 1,...,P
            % We express the correlation as a matrix multiplication
            % This matrix has a circulant structure (generated using gallery)
            % measurements(:,idx) = abs( gallery('circul', mask(:,idx)) * ...
            %                                                currentGuess );
            % The multiplication with a circulant matrix can be efficiently
            % implemented using FFTs
            % (discrete circular correlation/convolutions)         
            y(:,idMasks) = fft( fftFun .* maskFFT(:,idMasks) );
        end


        % Each pass (combination of HIO+ER) 
        for ipasses = 1:nPasses
        
        
        % In each HIO iteration iteration, ...
        for idx = 1:nHIO


            %% First projection - Set magnitudes using measurements
            phase = angle(y);
            P_1 = mag_meas .* exp(1i*phase);

            % Above is P_1(y) in Irene's terminology
            % Below is the first HIO projection
            w = 2*P_1 - y;
            

            %% Second projection - Project onto Range(A)
            % We may compute these effectively using FFTs

            % Pseudo-inverse operator
            projSoln = zeros(d,1);
            for idMasks = 1:nMasks
                projSoln = projSoln + ...
                    fft( (conj(maskFFT(:,idMasks))./normFac) .* ...
                                                ifft( w(:,idMasks) ) );

            end
            
            % Forward operator
            % For efficiency, compute these using FFTs
            % Recall that our measurements are correlations, which can be
            % efficiently evaluated using FFTs
            % Pre-store FFT of the current guess
            fftFun = ifft(projSoln);

            for idMasks = 1:nMasks
                % Measurements 
                % Structure of measurements - |corr(y,w)| = b^i, i = 1,...,P
                % We express the correlation as a matrix multiplication
                % This matrix has a circulant structure (generated using gallery)
                % measurements(:,idx) = abs( gallery('circul', mask(:,idx)) * ...
                %                                                currentGuess );
                % The multiplication with a circulant matrix can be efficiently
                % implemented using FFTs
                % (discrete circular correlation/convolutions)         
                P_RA(:,idMasks) = fft( fftFun .* maskFFT(:,idMasks) );
            end

            % This is P_2(w) in Irene's notation
            y_tilde = 2*P_RA - w;
            
            % Final update
            prevSoln = y;
            y = 0.5*(y + y_tilde);
            
            
            if( norm(prevSoln(:)-y(:)) <= pars.threshold )
                break;
            end
        end
        
        % Now, run iterations of ER
        % In each iteration, ...
        for idx = 1:nER


            %% First projection - Set magnitudes using measurements
            phase = angle(y);
            P_1 = mag_meas .* exp(1i*phase);

            % Above is P_1(y) in Irene's terminology
            % Below is the first HIO projection
            w = P_1;
            

            %% Second projection - Project onto Range(A)
            % We may compute these effectively using FFTs

            % Pseudo-inverse operator
            projSoln = zeros(d,1);
            for idMasks = 1:nMasks
                projSoln = projSoln + ...
                    fft( (conj(maskFFT(:,idMasks))./normFac) .* ...
                                                ifft( w(:,idMasks) ) );

            end
            
            % Forward operator
            % For efficiency, compute these using FFTs
            % Recall that our measurements are correlations, which can be
            % efficiently evaluated using FFTs
            % Pre-store FFT of the current guess
            fftFun = ifft(projSoln);

            for idMasks = 1:nMasks
                % Measurements 
                % Structure of measurements - |corr(y,w)| = b^i, i = 1,...,P
                % We express the correlation as a matrix multiplication
                % This matrix has a circulant structure (generated using gallery)
                % measurements(:,idx) = abs( gallery('circul', mask(:,idx)) * ...
                %                                                currentGuess );
                % The multiplication with a circulant matrix can be efficiently
                % implemented using FFTs
                % (discrete circular correlation/convolutions)         
                P_RA(:,idMasks) = fft( fftFun .* maskFFT(:,idMasks) );
            end

            % This is P_2(w) in Irene's notation
            y_tilde = P_RA;
            
            % Final update
            y = y_tilde;

        end
        
        end         % Done with a pass
        
        % Done - obtain image space solution
        % Pseudo-inverse operator
        projSoln = zeros(d,1);
        for idMasks = 1:nMasks
            projSoln = projSoln + ...
                fft( (conj(maskFFT(:,idMasks))./normFac) .* ...
                                            ifft( y(:,idMasks) ) );

        end


        
    case {'global-randGauss'}
        % These are global masks
        % Problem dimension
        d = pars.d;

        % We will work with absolute value (and not squared absolute value) of the
        % measurements
        measurements = sqrt(abs(measurements));
        
        % Max. iterations
        nIter = pars.altProjIters;

        % Start with initial guess
        projSoln = initGuess;

        % Pseudoinverse of measurement matrix
        pinv_measurementMat = pinv( maskEntries );

        % In each iteration, ...
        for idx = 1:nIter

            % Forward operator
            y = ( maskEntries * projSoln );

            % First projection - Set magnitudes using measurements
            phase = unwrap( angle(y) );
            y = measurements .* exp(1i*phase);

            % Second projection - Apply adjoint operator
            prevSoln = projSoln;                % Store previous iterate solution    
            projSoln = pinv_measurementMat * y;
            if( norm(prevSoln-projSoln) <= pars.threshold )
                return;
            end
        end

end

return
