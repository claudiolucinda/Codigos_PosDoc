function [mkup,cmg]=mg_costs(alpha,thetaNL,Delta,characteristic,Data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Baseline code for estimating Marginal Costs assuming Nash-Bertrand
% competition
% Claudio R. Lucinda
% 2017 - FEA-RP/USP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fnames=fieldnames(Data);
for i=1:length(fnames)
    eval([fnames{i} '=Data.' fnames{i} ';']);
end

[Ela_temp,~]=ElastNestedLogit(alpha,thetaNL,Delta,characteristic,Data);
meanval_start=deltaNL_PAR(thetaNL,Data);
[~,~,sj,~,~,~,~,~,~] = NLShareCalculation(thetaNL,meanval_start,Data);

mkup=zeros(size(sj,1),1);
cmg=zeros(size(sj,1),1);

for i=1:max(Data.cdid)
    marc=marca(Data.cdid==i,1);
    temp=unique(marc);
    dum_marca_t=zeros(size(marc,1),size(temp,1));
    for j=1:size(temp,1)
        dum_marca_t(marc==temp(j,1),j)=1;
    end
    own_mat=dum_marca_t*dum_marca_t';    
        
    
    El=Ela_temp(Data.cdid==i,Data.cdid==i)'.*own_mat;
    ss=sj(Data.cdid==i,1);
    tt=tax(Data.cdid==i,1);
    pp=price(Data.cdid==i,1);
    w=-El\ss;
    w=w./ss;
    mkup(Data.cdid==i,1)=w;
    cmg(Data.cdid==i,1)=(1-w).*(1-(tt./100)).*pp;
    
    
    
end    
