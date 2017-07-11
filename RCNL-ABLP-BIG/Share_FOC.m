function [sij,sijg,sj,sjg]=Share_FOC(oldprice,parameters,x1,x2,OData)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to return the shares and conditional shares 
% as a function of the prices and other stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arguments in OData
% gmmresid,nestid,dfull,vfull,income,pos

fnames=fieldnames(OData);
for i=1:length(fnames)
    eval([fnames{i} '=OData.' fnames{i} ';']);
end


nx1=size(x1,2);
nx2=size(x2,2);
theta2=parameters(nx1+1:nx1+nx2,:);
theta1=parameters(1:nx1,:);
Sigseg=parameters(nx1+nx2,:);

temp=unique(nestid);
temp_ooo=zeros(size(nestid,1),size(temp,1));
for j=1:size(temp,1)
    temp_ooo(nestid==temp(j,1),j)=1;
end

%temp_ooo=sparse(dummyvar(nestid));

% Assembling stuff for using back mufunc

x1(:,pos)=oldprice./income;
Pack.x1=x1;
Pack.x2=x2;
Pack.dfull=dfull;
Pack.vfull=vfull;
%ns=size(dfull,2);
Pack.ns=ns;

little_mu=mufunc(theta2,Pack);
little_delta=x1*theta1+gmmresid;

ind_util=little_mu+(little_delta)*ones(1,size(dfull,2));
numer1=exp((ind_util)./(1-Sigseg));

Big_selector=temp_ooo*temp_ooo';
denom1=Big_selector*(numer1);

numer2=denom1.^(1-Sigseg);

denom2=zeros(size(denom1));

for j=1:ns
    temp=sum(unique(numer2(:,j)),1);
    denom2(:,j)=temp*ones(size(denom2(:,j)));
end

denom2=1+denom2;

sijg=numer1./denom1;
sijg(isnan(sijg)) = 0;

sig=numer2./denom2;
sig(isnan(sig))   = 0;

sij=sig.*sijg;
sij(isnan(sij)) = 0;

sjg=mean(sijg,2);
sj=mean(sij,2);








