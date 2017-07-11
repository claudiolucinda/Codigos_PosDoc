function [foc]=Foc_SYS(oldprice,parameters,x1,x2,OData)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that returns the FOC for a Bertrand equilibrium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arguments in OData
% gmmresid,nestid,dfull,vfull,income,pos,tax,mg_cost

fnames=fieldnames(OData);
for i=1:length(fnames)
    eval([fnames{i} '=OData.' fnames{i} ';']);
end

[~,~,spoint,~] = Share_FOC(oldprice,parameters,x1,x2,OData);
elapoint       = Elast_FOC(oldprice,parameters,x1,x2,OData);

mkup=((1-tax).*oldprice-mg_cost)./((1-tax).*oldprice);
w=spoint.*mkup;
foc=spoint+elapoint*w;