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

function [xrec, etime] = postprocess_HIOER( x, measurements, m, d, K, L, varargin )


%% Initialization

% process optional arguments
if( nargin>6 )
    nits = varargin{1};
    nIter = nits(1);
    nHIO = nits(2); nER = nits(3);
else
    % Default no. of HIO vs ER iterations
    nHIO = 25; nER = 5;        % 20 iter. of HIO + 10 iter. of ER     
    nIter = 60;
end

if (nargin>7 )
    M = varargin{2};
    Mproj = varargin{3};
end

tic;

% We will work with absolute value (and not squared absolute value) of the
% measurements
mag_meas = sqrt(measurements(:));

% threshold
threshold = 1e-14;

% Store result of applying forward operator/measurement matrix here
% y = zeros(K, d);
% P_RA = zeros(K, d);

% % Normalization factor for aplpying the (pseudo) inverse transform
% normFac = sum( abs(maskFFT).^2, 2 );

% Max. iterations
nPasses = round(nIter/(nHIO+nER));

% Normalization for pseudoinverse
% normFac = zeros(d,1);
% for l = 0:K-1
%     shft = l*d/K;
%     normFac = normFac + circshift(abs(m).^2,shft);
% end
% normFac = normFac(1);

% for l = 0:d-1
%     normFac = normFac + circshift(abs(m).^2, l);
% end

% % Sub-sampled Fourier matrix
% Fd = dftmtx(d);
% Fk = Fd(1:d/K:end,:);
% clear Fd;
% 
% 
% % Measurement matrix
% M = zeros(L*K, d);
% ctr = 1;
% for ix = 1:d/L:d
%     M((ctr-1)*K+1:ctr*K, :) = Fk*diag(circshift(m,ix-1));
%     ctr = ctr+1;
% end
% Mproj = M*pinv(M);

% Initial guess
x = full(x);

%% Forward operator
% for l=0:d/L:d-1
%     tmp = fft(x.*circshift(m,l), d);
%     y(:,l*(d/L)+1) = tmp(1:d/K:end);
% end

y = M*x; 

% for l = 1:L
%     % Measurements 
%     shift = (l-1)*d/L;
%     y(:,l) = fft(x.*circshift(m,shift), d);
% end

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
%     xrec = zeros(d,1);
%     for l = 1:L
%         shift = (l-1)*d/L;
%         xrec = xrec + circshift(conj(m),shift) .* ifft( w(:,l) ) ./ normFac;
%     end
%     for l = 1:d
%         xrec = xrec + circshift(conj(m), l-1) .* ...
%                                     ( Fk'*w(:,l) ) ./ normFac;
%     end
%     
%     % Forward operator
%     for l=0:d-1
%         tmp = fft(xrec.*circshift(m,l), d);
%         P_RA(:,l+1) = tmp(1:d/K:end);
%     end
    
%     P_RA = reshape( Mproj*w(:), K, d);
    P_RA = Mproj*w;

%     for l = 1:L
%         shift = (l-1)*d/L;
%         P_RA(:,l) = fft(xrec.*circshift(m,shift), d);
%     end

    % This is P_2(w) in Irene's notation
    y_tilde = 2*P_RA - w;

    % Final update
    prevSoln = y;
    y = 0.5*(y + y_tilde);


    if( norm(prevSoln(:)-y(:)) <= threshold )
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

%     % Pseudo-inverse operator
%     xrec = zeros(d,1);
% %     for l = 1:L
% %         shift = (l-1)*d/L;        
% %         xrec = xrec + circshift(conj(m),shift) .* ifft( w(:,l) ) ./ normFac;
% %     end
%     for l = 1:d
%         xrec = xrec + circshift(conj(m), l-1) .* ...
%                                     ( Fk'*w(:,l) ) ./ normFac;
%     end
%     
%     % Forward operator
%     for l=0:d-1
%         tmp = fft(xrec.*circshift(m,l), d);
%         P_RA(:,l+1) = tmp(1:d/K:end);
%     end
    
%    P_RA = reshape( Mproj*w(:), K, d);
    P_RA = Mproj*w;
    
%     for l = 1:L
%         shift = (l-1)*d/L;
%         P_RA(:,l) = fft(xrec.*circshift(m,shift), d);
%     end

    % This is P_2(w) in Irene's notation
    y_tilde = P_RA;

    % Final update
    y = y_tilde;

end

end         % Done with a pass

% Done - obtain image space solution
% % Pseudo-inverse operator
% xrec = zeros(d,1);
%     for l = 1:d
%         xrec = xrec + circshift(conj(m), l-1) .* ...
%                                     ( Fk'*w(:,l) ) ./ normFac;
%     end
% % for l = 1:L
% %     shift = (l-1)*d/L;    
% %     xrec = xrec + circshift(conj(m),shift) .* ifft( w(:,l) ) ./ normFac;
% % end

    xrec = inv(M'*M)*M'*w(:);

etime = toc;

return
