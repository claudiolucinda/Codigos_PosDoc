function [sij sijg sj sjg s0 numer1 denom1 numer2 denom2]=NLShareCalculation(theta,Delta,Data)

% NLShareCalculation - Nested Logit share calculation
% Originally from Grigolon and Verboven 2014 REStat
% Adapted by Claudio Lucinda

% Syntax:   [sij sijg sj sjg s0 numer1 denom1 numer2 denom2]=NLShareCalculation(thetaNL,delta,Data)
%
% Inputs:
%    input1 - Non linear parameters - Including rho
%    input2 - Vector of Mean Valuations
%    input2 - Data
%
% Outputs:
%    sij    = individual market shares
%    sijg   = individual conditional share of choosing good j from group g
%    sj     = market shares
%    sjg    = conditional share of choosing good j from group g 
%    s0     = outside good
%    numer1 = numerator1
%    denom1 = denominator1
%    numer2 = numerator 2
%    denom2 = denominator 2

% Author: Laura Grigolon and Frank Verboven
% August 2012;

% Unpack

warning('off','MATLAB:nearlySingularMatrix');

%global args

fnames=fieldnames(Data);
for i=1:length(fnames)
    eval([fnames{i} '=Data.' fnames{i} ';']);
end


% prods           = Data.prods;
% nmkt            = Data.nmkt;
% Nnodes          = Data.Nnodes;
% qweightrprods   = Data.qweightrprods;
% multidumseg     = Data.multidumseg;
% xvu             = Data.xvu;

mu = mufunc(theta,Data);
Sigseg = theta(end);
expDelta = Delta*ones(1,ns);

% Market share calculation
%mudel           = (bsxfun(@plus,delta,xvu.*RandCoeff)) ./(1-Sigseg);
numer1          = exp((mu+expDelta)./(1-Sigseg));

temp_oo=sparse(dummyvar(cdid));
temp_ooo=sparse(dummyvar(nestid));
matrix_sel1=sparse(temp_oo*temp_oo');

tam=0;
for i=1:max(Data.cdid)
    tam=tam+(size(Data.cdid(Data.cdid==i,:),1)).^2;
end

matrix_sel2=spalloc(size(matrix_sel1,1),size(matrix_sel1,1),tam);


for i=1:max(Data.cdid)
   ttemp_ooo=temp_ooo(Data.cdid==i,:); 
   matrix_sel2(Data.cdid==i,Data.cdid==i)=ttemp_ooo*ttemp_ooo'; 
    
end    
%matrix_sel2=sparse(temp_ooo*temp_ooo');
Big_selector=matrix_sel1.*matrix_sel2;
denom1=Big_selector*(numer1);

numer2=denom1.^(1-Sigseg);

denom2=zeros(size(denom1));

% Versão sem paralelização
for i=1:max(Data.cdid)
    %Sel_denom2=denom1(cdid==i,:);
    denom_temp2=zeros(size(denom1(Data.cdid==i,:)));
    for j=1:ns
        temp=sum(unique(numer2(Data.cdid==i,j)),1);
        denom_temp2(:,j)=temp*ones(size(denom_temp2(:,j)));
    end
    denom2(Data.cdid==i,:)=denom_temp2;
end

denom2=1+denom2;

sijg=numer1./denom1;
sijg(isnan(sijg)) = 0;

sig=numer2./denom2;
sig(isnan(sig))   = 0;

sij=sig.*sijg;
sij(isnan(sij)) = 0;

sjg=mean(sijg,2);
sg=mean(sig,2);
sj=mean(sij,2);

s0=mean(1-matrix_sel1*sij,2);

% Versão com paralelização
% 
% 
% parfor i=1:max(Data.cdid)
%     %ctemp=cdid;
%     Sel_I_ijg=I_ijg(Data.cdid==i,:);
%     I_temp2=zeros(size(Sel_I_ijg));
%     for j=1:Data.ns
%         I_temp=sum(unique(Sel_I_ijg(:,j)),1);
%         I_temp2(:,j)=exp(I_temp);
%     end
%     I_i2{i}=I_temp2;
% end
% 
% 
% I_i2=real(cat(1,I_i2{:}));



%D_sel2=dummyvar([cdid nestid]);

%temp=pinv(D_sel2);
%Proj_mat=D_sel2*((D_sel2'*D_sel2)\D_sel2);

%I_i=Proj_mat*I_ijg;
% I_i=log(1+I_i2);
% 
% 
% sijg=exp(numer1)./denom1;
% 
% sij = ((numer1.*exp(I_ijg))./(exp(I_ijg./(1-Sigseg)))).*(exp(I_ijg)./exp(I_i));
% sij(isnan(sij)) = 0;
% 
% sg = Big_selector*sij;
% 
% sijg = sij./sg;
% sijg(isnan(sijg)) = 0;
% 
% sj=mean(sij,2);
% sj(isnan(sj)) = 0;
% 
% 
% sjg=mean(sijg,2);
% sjg(isnan(sjg)) = 0;
% 
% s0=sum(1-matrix_sel1*sij,2);

end