function se     = seNL(parameters,Data)

nlin  = Data.nlin;
nrcNL = Data.nrcNL;

theta2   = parameters(nlin+1:nlin+nrcNL,:);
deltaopt = deltaNL(theta2,Data);
derdel   = jacobNL(theta2,deltaopt,Data);
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