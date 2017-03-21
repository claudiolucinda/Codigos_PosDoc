function [f df]=gmmNL(thetaNL,BLPdata)
% gmmNL - This function is used to calculate the GMM estimator for the
% Nested Logit function
% 
% [f df]=gmmNL(thetaNL,BLPdata)
%
% Inputs:
%    input1 - Non linear parameters nested logit
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
d   = deltaNL(thetaNL,BLPdata);
%% GMM
if max(isnan(d)) == 1
 f  = 1e10;	  
else 
% Use relationship between linear and non-linear parameters from step 4,
% resulting from the FOC's
    betaNL  = BLPdata.invxzwzx*(BLPdata.xzwz*d);
    % error term
    csi     = d - BLPdata.Xexo*betaNL;
	f       = csi'*BLPdata.Z*BLPdata.W*BLPdata.Z'*csi;
    save betaNL betaNL;
end
%% Gradient
if max(isnan(d)) == 1
    df=1e10;
elseif nargout>1
        temp    = jacobNL(thetaNL,d,BLPdata);
        df      = 2*temp'*BLPdata.Z*BLPdata.W*BLPdata.Z'*csi;
end

end

