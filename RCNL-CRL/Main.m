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
cd  'C:\Users\claudiolucinda\Dropbox\GIT-Posdoc\RCNL-CRL\'
data_dir='C:\Users\claudiolucinda\Dropbox\Pós Doc 2016\Paper\Data\';

% IMportando os dados
regressors=dlmread([data_dir 'regressors_T.txt'],'\t',1,0);
other=dlmread([data_dir 'extra_data_T.txt'],'\t',1,0);
other2=dlmread([data_dir 'extra_data2_T.txt'],'\t',1,0);

path(path,'C:/Users/claudiolucinda/Dropbox/GIT-Posdoc/RCNL-CRL/');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IV Regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mean_inv=other(:,end-1);
cond_sh=other(:,end);
regrs=[ones(size(cond_sh)) regressors(:,2:end) cond_sh];
z=other(:,1:6);
incl_inst=[ones(size(cond_sh)) regressors(:,2:end-1)];
IV=[incl_inst z];
beta_IV_ini=ivregression(mean_inv,regrs,IV);

rho_IV_ini=beta_IV_ini(end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting up stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x2=[regressors(:,2) regressors(:,3) regressors(:,4) cond_sh]; 
x1=[ones(size(cond_sh)) regressors(:,3:end-1)];
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
theta2_test=[beta_IV_ini(2)./2.3043e+05;1e-5;1e-5;rho_IV_ini];

% Testing MUFUNC
% test_mu=mufunc(theta2_test,Data);
% Testing OK
Delta_Init=regrs(:,1:end-1)*beta_IV_ini(1:end-1);

% Testing NLShareCalculation
Data.cdid=cdid;
Data.nestid=nestid;

% tic
% [sij,bb,sh,dd,ee] = NLShareCalculation(theta2_test,Delta_Init,Data);
% toc

% Testing OK
% Tempo (64bit architecture, large cities): 19.4136 seconds
% Tempo paralelizado: 24.879 seconds ???

% Testing deltaNL
delta0NL=Delta_Init;
save([data_dir 'mvaloldNL.mat'],'delta0NL');
Data.data_dir=data_dir;
Data.sj=s_jt;
% tic
% jjj=deltaNL(theta2_test,Data);
% toc

% Número de Iterações: 13
% Elapsed time is 250.093450 seconds.

tic
jjj=deltaNL_PAR(theta2_test,Data);
toc

% Um por vez: Elapsed time is 45.901219 seconds. 
% Paralelização: Elapsed time is 20.186337 seconds.
% Testing OK
Data.IV=IV;
Data.x1=x1;

% 
tic
kkk=gmmNL(theta2_test,Data);
toc

tic
kkk=jacobNL(delta0NL,theta2_test,Data);
toc
