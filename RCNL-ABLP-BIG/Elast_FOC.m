function [elast]=Elast_FOC(oldprice,parameters,x1,x2,OData)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to compute the price elasticity of demand as a function of prices
% and a whole bunch of stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arguments in OData
% gmmresid,nestid,dfull,vfull,income,pos

fnames=fieldnames(OData);
for i=1:length(fnames)
    eval([fnames{i} '=OData.' fnames{i} ';']);
end

temp=unique(nestid);
temp_ooo=zeros(size(nestid,1),size(temp,1));
for j=1:size(temp,1)
    temp_ooo(nestid==temp(j,1),j)=1;
end


nx1=size(x1,2);
nx2=size(x2,2);
theta2=parameters(nx1+1:nx1+nx2,:);
theta1=parameters(1:nx1,:);
Sigseg=parameters(nx1+nx2,:);

[sij,sijg,sj,~]=Share_FOC(oldprice,parameters,x1,x2,OData);

priceinc=oldprice./income;

part3           = (sij*sij');
part2           = ((Sigseg./(1-Sigseg)).*sijg) * sij';
part1           = diag((sum(sij,2)./(1-Sigseg)));
%temp_ooo        = dummyvar(nestid);
Big_selector    = temp_ooo*temp_ooo';
dsD = (part1 - Big_selector.*part2 - part3)./ns;
PriceOverShare  = (priceinc./sj)*ones(1,size(dsD,2));
elast = alpha2.*PriceOverShare.*dsD;



