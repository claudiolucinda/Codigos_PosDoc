%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN FILE
% RC-NL
% Claudio R. Lucinda
% FEA-RP/USP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the number of cores
ncores=str2double(getenv('NUMBER_OF_PROCESSORS'))-1;
matlabpool close force local
matlabpool('open',ncores)

clear; clc;
rng(89282250);

warning('off','MATLAB:nearlySingularMatrix');
cd  'C:\Users\claudiolucinda\Dropbox\GIT-Posdoc\RCNL-ABLP-BIG\'
data_dir='C:\Users\claudiolucinda\Dropbox\Pós Doc 2016\Paper\Data\';

% IMportando os dados
regressors=dlmread([data_dir 'regressors_T_M.txt'],'\t',1,0);
regressors(:,2)=regressors(:,2)./10;
regressors(:,3)=regressors(:,3)./1000;

other=dlmread([data_dir 'extra_data_T_M.txt'],'\t',1,0);
other2=dlmread([data_dir 'extra_data2_T_M.txt'],'\t',1,0);

cdid=other(:,7);
nmkts=max(cdid);

dum_mkt=spalloc(size(cdid,1),nmkts,size(cdid,1));
for i=1:nmkts
    dum_mkt(cdid==i,i)=1;
end

dum_mkt(:,1)=[];

marca=other2(:,end-1)+1;

dum_marca=sparse(dummyvar(marca));
%sel=full(sum(dum_marca,1))<10000;
%dum_marca(:,sel)=[];
dum_marca(:,[1 8 9 18 26 27])=[];

pmat_b=[dum_mkt dum_marca];
%pmat_b=[dum_marca(:,2:end)];
%pmat_b=[dum_mkt];

path(path,'C:/Users/claudiolucinda/Dropbox/GIT-Posdoc/RCNL-ABLP-BIG/');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IV Regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mean_inv=other(:,end-1);
cond_sh=other(:,end);
regrs=sparse([ones(size(cond_sh)) regressors(:,2:end) pmat_b cond_sh]);
z=other(:,1:6);
incl_inst=sparse([ones(size(cond_sh)) regressors(:,2:end-1) pmat_b]);
IV=[incl_inst z];
beta_IV_ini=ivregression(mean_inv,regrs,IV);

norm_RC_ini=abs(.5.*beta_IV_ini(3:4));
%norm_RC_ini=abs(.5.*beta_IV_ini(4));
rho_IV_ini=beta_IV_ini(end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting up stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x2=sparse([regressors(:,2) regressors(:,3) regressors(:,4) cond_sh]); 
%x2=sparse([regressors(:,2) regressors(:,4) cond_sh]); 
x1=sparse([ones(size(cond_sh)) regressors(:,3:end) pmat_b]);
s_jt=regressors(:,1);
ns=250;
nobs=size(x1,1);
ncoefs=size(x1,2)+size(x2,2);
K = size(x2,2);
cdid=other(:,7);
nestid=other2(:,end)+1;
cdindex=find(diff(cdid));
cdindex=[cdindex;size(cdid,1)];
nmkts=max(cdid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating the random draws
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stored_draws=0;
v=[];
if stored_draws==1
    v=load([data_dir 'sigdraws2.mat']);
    v=v.v;
else
    for k=1:K
        vtmp=random('Normal',0,1,max(nmkts),ns/2);
        v=[v vtmp -vtmp];
    end
    save([data_dir 'sigdraws2.mat'],'v');
    clear vtmp
end
demogr=dist_kms(ns,nmkts);
vfull=v(cdid,:);
dfull=demogr(cdid,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data.dfull=dfull;
Data.vfull=vfull;
Data.ns=ns;
Data.x2=x2;
nnorm_RC_ini=beta_IV_ini(2)./mean(dfull(1,:),2);
theta2_test=[nnorm_RC_ini;norm_RC_ini;rho_IV_ini];

% Testing MUFUNC
%test_mu=mufunc(theta2_test,Data);
% Testing OK
Delta_Init=[ones(size(cond_sh)) regrs(:,3:end-1)]*[beta_IV_ini(1);beta_IV_ini(3:end-1)];

% Testing NLShareCalculation
Data.cdid=cdid;
Data.nestid=nestid;


 
% 
% % Testing OK
% % Tempo (64bit architecture, large cities): 19.4136 seconds
% % Tempo paralelizado: 24.879 seconds ???
% 
% % Testing deltaNL
delta0NL=Delta_Init;
save([data_dir 'mvaloldNL.mat'],'delta0NL');
Data.data_dir=data_dir;
Data.sj=s_jt;
Data.lshare=log(s_jt);
meanval_start=deltaNL_PAR(theta2_test,Data);
%[sij,sijg,sj,sjg,s0,numer1,denom1,numer2,denom2] = NLShareCalculation(theta2_test,meanval_start,Data);
%teste=[s_jt(cdid==1,:) sj(cdid==1,:)];


Data.IV=IV;
Data.x1=x1;
W = (IV'*IV) \ eye(size(IV,2));
Data.W=W;
%f=gmmNL(theta2_test,Data);
%[f,df]=gmmNL(theta2_test,Data);

%%

options = optimset( 'Display','iter',...
    'GradObj','on','TolCon',1E-6,...
    'TolFun',1E-6,'TolX',1E-6,...
    'Hessian', 'off','DerivativeCheck','off','FinDiffType','central');

%     [thetaRCNL, FctValRCNL, exitflagRCNL, output, gradient] = ...
%         fminunc(angmmNL,theta20RCNL,options);

%
diary off
obj_fun=@(theta2)gmmNL(theta2,Data);

% x_L     = [-1e3;0;0];           % lower bound (only for the nested logit, to avoid numerical problems when pho gets negative)
% x_U     = [1e3;4e2;0.975];       % upper bound (only for the nested logit, to avoid numerical problems when pho gets close to 1)

x_L     = [1.5*theta2_test(1);0;0;0];           % lower bound (only for the nested logit, to avoid numerical problems when pho gets negative)
x_U     = [1.5*abs(theta2_test(1:3));0.975];       % upper bound (only for the nested logit, to avoid numerical problems when pho gets close to 1)
[thetaRCNL, FctValRCNL, exitflagRCNL,~,~,~,hessRCNL] = ...
    knitromatlab(obj_fun,theta2_test, [], [], [], [], x_L, x_U, [], [],options, 'knitroBLP.opt');

% First Stage Results
meanval_FINAL=deltaNL_PAR(thetaRCNL,Data);
theta1=(x1'*IV * W * IV'*x1)\(x1'*IV * W * IV'*meanval_FINAL);
gmmresid = meanval_FINAL - x1*theta1;

save([data_dir 'gmmresid.mat'],'gmmresid');
save([data_dir 'mvaloldNL.mat'],'meanval_FINAL');

% Second Stage Results
[bhat,uhat]=ivregression(meanval_FINAL,x1,IV,W);

g=bsxfun(@times,uhat,IV);
gstar=bsxfun(@minus,g,mean(g))./size(x1,2);
S=gstar'*gstar;

delta0NL=meanval_FINAL;
save([data_dir 'mvaloldNL.mat'],'delta0NL');

%W2=(size(x1,2).^2*S)\eye(size(IV,2));
%Data.W=W2;

th12RCNL2S      = [theta1 ; thetaRCNL];
sterrRCNL     = seNL(th12RCNL2S,Data);

% Random Coefficient Nested Logit wrt X1
alpha         = theta1(17);
priceinc      = x1(:,17);

[ElaRCNL,DiversionRCNL] = ElastNestedLogit(alpha,thetaRCNL,delta0NL,priceinc,Data);

meanElast=mean(diag(ElaRCNL),1);

