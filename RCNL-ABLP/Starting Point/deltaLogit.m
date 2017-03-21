function [delta1, i] = deltaLogit(theta,Data)

% DeltaNL - Contraction Mapping for Logit
% Syntax:  f = deltaLogit(theta,Data)
%
% Inputs:
%    input1 - Non linear parameters
%    input2 - Data
%
% Outputs:
%    Mean Value Delta
%
% MAT-files required: mvalold
% Subfunctions: LogitShareCalculation
%

% Author: Laura Grigolon and Frank Verboven
% August 2012;

load mvalold
k   = 100;
km  = 1e-14;
i   = 0;
% if max(isnan(delta0))  == 1 || isreal(delta0) == 0
if max(isnan(delta0))  == 1
    deltastart         =  zeros(Data.nobs,1);
else
    deltastart         = delta0;
end

% Unpack
logobsshare            = Data.logobsshare;

while k > km
    % Market Share
    if numel(theta) == 1
        [sh, ~] = LogitShareCalculation(theta,deltastart,Data);
    elseif numel(theta) == 2
        [sh, ~] = LogitShareCalculation2RC(theta,deltastart,Data);
    end
    sh      = reshape(sh,[Data.nobs 1]);   % 2D market shares
    % Contraction Mapping
    delta1  = deltastart + logobsshare-log(sh);
    if max(isnan(delta1)) == 1
%         disp('No Convergence - delta calculation failed: OVERFLOW')
    break
    end
    
    if isreal(delta1) == 0
        break
    end
    
    i = i + 1;
    if i > 2500
%         disp('No Convergence - delta convergence failed')
        break
    else
    k           = max(abs(delta1-deltastart));
    deltastart  = delta1;
    end
end

% disp(['# of iterations for delta convergence:  ' num2str(i)])

if isreal(delta1)   == 0
    delta0          =  zeros(Data.nobs,1);
else
    delta0 = delta1;
end

save mvalold delta0