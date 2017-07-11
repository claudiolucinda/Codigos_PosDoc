function [newprice, newshare, checker]=Solver_BERT(oldprice,parameters,x1,x2,Data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that returns the new price vector and new share vector 
% As a Nash Bertrand Equilibrium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fnames=fieldnames(Data);
for i=1:length(fnames)
    eval([fnames{i} '=Data.' fnames{i} ';']);
end

newprice=zeros(size(oldprice));
newshare=zeros(size(oldprice));
checker=zeros(size(unique(cdid)));

% Cenário 1
% Versão sem paralelização
% solving mkt by mkt
options=optimset('UseParallel','always','Algorithm','trust-region-reflective','Display','final','MaxFunEvals',1e4,'MaxIter',1e2,'TolX',1e-6,'TolFun',1e-6);

for i=1:max(Data.cdid)
    selector=(mg_cost>.01 & cdid==i);
    OData.mg_cost=mg_cost(selector,:);
    OData.gmmresid=gmmresid(selector,:);
    OData.tax=tax(selector,:)./100;
    OData.pos=pos;
    OData.dfull=dfull(selector,:);
    OData.vfull=vfull(selector,:);
    OData.nestid=nestid(selector,:);
    OData.income=income(selector,:);
    OData.ns=ns;
    OData.alpha2=alpha2;
    startprice=oldprice(selector,:);
    x1mkt=x1(selector,:);
    x2mkt=x2(selector,:);
    disp(['Solving Market: ' num2str(i)]);

    [pfin,~,~,exflag]=lsqnonlin(@(oldprice) Foc_SYS(oldprice,parameters,x1mkt,x2mkt,OData),startprice,zeros(size(startprice)),[],options);
    [~,~,sfin,~]=Share_FOC(pfin,parameters,x1mkt,x2mkt,OData);
    newprice(selector,:)=pfin;
    newshare(selector,:)=sfin;
    checker(i,:)=exflag;
end

