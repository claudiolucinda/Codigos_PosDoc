function se=seRCLogit(parameters,Data)
% se - Standard Errror calculation of the RC Logit
% Syntax:   se=seRCLogit(parameters,Data)
%
% Inputs:
%    input1 - Non linear parameters
%    input2 - Data
%
% Outputs:
%    Standard Errors
%
%
% Subfunctions: deltaLogit, jacobLogit
%

% Author: Laura Grigolon
% August 2012;

nlin = Data.nlin;
nrc  = Data.nrc;

theta2   = parameters(nlin+1:nlin+nrc,1);
deltaopt = deltaLogit(theta2,Data);
derdel   = jacobLogit(theta2,deltaopt,Data);
derksi   = [-Data.Xexo derdel]'*Data.Z;
vv       = (derksi*Data.W*derksi')\eye(size(derksi,1));
xi       = deltaopt-Data.Xexo*parameters(1:nlin,:);

covg     = zeros(size(Data.Z,2));
for ii   = 1:length(Data.Z),
    covg = covg + Data.Z(ii,:)'*Data.Z(ii,:)*(xi(ii)^2);
end

varcovar = vv*derksi*Data.W*covg*Data.W*derksi'*vv;
se       = sqrt(diag(varcovar));
end