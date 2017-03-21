function [f, df]     = gmmLogit(theta,BLPdata)
% gmmLogit - This function is used to calculate the GMM estimator for the
% Logit function
% 
% Syntax: [f df]     = gmmLogit(theta,BLPdata)
%
% Inputs:
%    input1 - Non linear parameters
%    input2 - Data
%
% Outputs:
%    f = gmm function value; df = gradient of f 
%
% Subfunctions: deltaLogit(theta,BLPdata) ; jacobLogit
%

% Author: Laura Grigolon and Frank Verboven
% August 2012;

%% Contraction Mapping
[d,iterat]   = deltaLogit(theta,BLPdata);

%% GMM
if max(isnan(d)) == 1 || iterat>2500 || isreal(d)==0
 f  = 1e10;
else 
 
% Use relationship between linear and non-linear parameters from step 4,
% resulting from the FOC's
    beta = BLPdata.invxzwzx*(BLPdata.xzwz*d);
    % error term
    csi  = d - BLPdata.Xexo*beta;
	f    = csi'*BLPdata.Z*BLPdata.W*BLPdata.Z'*csi;
    save beta beta;
end
%% Gradient
if max(isnan(d)) == 1 || iterat>2500 || isreal(d)==0    
    % isreal(d)==0 can happen because of the negative weights - this is a problem of quadrature methods
    df=1e10;
else
    temp    = jacobLogit(theta,d,BLPdata);
    df      = 2*temp'*BLPdata.Z*BLPdata.W*BLPdata.Z'*csi;
end

end

