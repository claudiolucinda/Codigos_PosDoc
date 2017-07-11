function [lsum]=Logsum_RCNL(oldprice, parameters, x1, x2, alpha2, OData)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that returns the logsum by market
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fnames=fieldnames(OData);
for i=1:length(fnames)
    eval([fnames{i} '=OData.' fnames{i} ';']);
end

nx1=size(x1,2);
nx2=size(x2,2);
theta2=parameters(nx1+1:nx1+nx2,:);
theta1=parameters(1:nx1,:);
Sigseg=parameters(nx1+nx2,:);

x1(:,pos)=oldprice./income;
lsum=zeros(size(unique(cdid)));

for i=1:max(cdid)
    selector=(mg_cost>.01 & cdid==i);
    gmmresidmkt=gmmresid(selector,:);
    Extradat.dfull=dfull(selector,:);
    Extradat.vfull=vfull(selector,:);
    Extradat.ns=ns;
    nestidmkt=nestid(selector,:);
    temp=unique(nestidmkt);
    temp_ooo=zeros(size(nestidmkt,1),size(temp,1));
    for j=1:size(temp,1)
        temp_ooo(nestidmkt==temp(j,1),j)=1;
    end

    x1mkt=x1(selector,:);
    Extradat.x2=x2(selector,:);
    little_mu=mufunc(theta2,Extradat);
    little_delta=x1mkt*theta1+gmmresidmkt;
    ind_util=little_mu+(little_delta)*ones(1,size(dfull,2));
    numer1=exp((ind_util)./(1-Sigseg));
    Big_selector=temp_ooo*temp_ooo';
    denom1=(1-Sigseg).*log(Big_selector*(numer1));
    ls=zeros(1,ns);
    incmkt=income(selector,:);
    for j=1:ns
        temp=sum(exp(unique(denom1(:,j))),1);
        ls(:,j)=log(temp);
        
    end
    lsum(i,:)=mean(ls,2)./(abs(alpha2)./mean(incmkt,1));





    
    
end