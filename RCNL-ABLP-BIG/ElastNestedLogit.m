function [ElaNLogit Diversion]=ElastNestedLogit(alpha,thetaNL,Delta,characteristic,Data)

% ElastNestedLogit - Calculate price elasticities for Nested Logit models
% Syntax:  f=ElastLogit(alpha,theta,delta,Data)
%
% Inputs:
%    input1 - alpha     = price coefficient
%    input2 - thetaNL   = random coefficient
%    input1 - delta     = mean value
%    input2 - Data
%
% Outputs:
%    prods x prods elasticity matrix
%
% Subfunctions: NLShareCalculation 
%
% Author: Laura Grigolon and Frank Verboven
% August 2012;
fnames=fieldnames(Data);
for i=1:length(fnames)
    eval([fnames{i} '=Data.' fnames{i} ';']);
end

Sigseg =thetaNL(end);

K=size(thetaNL,1);

[sij,sijg,sj,~,~,~,~,~,~] = NLShareCalculation(thetaNL,Delta,Data);


% Computing the derivative of share w.r.t. sigma

%mu=mufunc(thetaNL,Data);
tam=0;
for i=1:max(Data.cdid)
    tam=tam+(size(Data.cdid(Data.cdid==i,:),1)).^2;
end

derShareDelt=spalloc(size(cdid,1),size(cdid,1),tam);
ElaNLogit=spalloc(size(cdid,1),size(cdid,1),tam);
Diversion=spalloc(size(cdid,1),size(cdid,1),tam);

for i=1:max(Data.cdid)
    % derivative of share w.r.t. DELTA
    part3           = (sij(Data.cdid==i,:)*sij(Data.cdid==i,:)');
    part2           = ((Sigseg./(1-Sigseg)).*sijg(Data.cdid==i,:)) * sij(Data.cdid==i,:)';
    part1           = diag((sum(sij(Data.cdid==i,:),2)./(1-Sigseg)));
    temp_ooo        = sparse(dummyvar(Data.nestid(Data.cdid==i,:)));
    Big_selector    = temp_ooo*temp_ooo';
    dsD = (part1 - Big_selector.*part2 - part3)./Data.ns;
    derShareDelt(Data.cdid==i,Data.cdid==i)    = dsD;
    PriceOverShare  = (mean(characteristic(Data.cdid==i,:)./sj(Data.cdid==i,:),2))*ones(1,size(dsD,2));
    el = alpha.*PriceOverShare.*dsD;
    ElaNLogit(Data.cdid==i,Data.cdid==i)       = el;
    div = dsD.*(1./(diag(dsD)*ones(1,size(dsD,1))));
    Diversion(Data.cdid==i,Data.cdid==i)       = div;
    
    
end

end
