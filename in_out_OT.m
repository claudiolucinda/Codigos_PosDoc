    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arquivo Input-Output
% Copyright 2016 Cláudio R. Lucinda
% FEA-RP/USP
% Versão DEZEMBRO DE 2016
% Agora puxando diretamente do Pharmaeco os dados
% Cenário Brasil vs. GCR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(RandStream('mcg16807', 'Seed',0));
% Se usando processamento paralelo descomentar aqui
matlabpool close force local
matlabpool open 2

clear; clc;
rng(89282250);

cd  'C:\Users\claudiolucinda\Dropbox\GIT-Posdoc\'
data_dir='C:\Users\claudiolucinda\Dropbox\Pós Doc 2016\Paper\Data\';

% IMportando os dados
regressors=dlmread([data_dir 'regressors.txt'],'\t',1,0);
other=dlmread([data_dir 'extra_data.txt'],'\t',1,0);





ns=250;

y=other(:,end);
regrs=regressors(:,2:end);
x=[ones(size(y)) regrs];
x2=[regressors(:,3) regressors(:,4)]; %x(:,2)  ones(size(x,1),1)];
x1=x(:,1:end-1);
s_jt=regressors(:,1);

%x2(:,1)=x2(:,1)/1000;
%x1(:,end)=x1(:,end)/1000;

nobs=size(x,1);
ncoefs=size(x,2)+size(x2,2);

[n,K] = size(x2);

cdid=other(:,7);
cdindex=find(diff(cdid));
cdindex=[cdindex;size(cdid,1)];
nmkts=max(cdid);


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
demogr=[];
vfull=v(cdid,:);
dfull=[];
%dfull=demogr(cdid,:);
argums.dfull=dfull;

% IV Stuff Just to check 
z=other(:,1:7);
IV=[x(:,1:end-1) z];
% Starting Values - Variáveis não RC
miolo=pinv(IV);
tmp1=x'*IV;
t=(tmp1*miolo*x)\(tmp1*miolo*y);
% KindaOff, mas pq aqui não tô fazendo o robust


z=other(:,1:7);



IV=[x1(:,1:end-1) z];


% Starting Values variável RC - preço

theta2w=[1e-5 1e-5]'; % eps eps]';

[theti, thetj, theta2]=find(theta2w);

% Starting Values - Variáveis não RC
miolo=pinv(IV);

%invA = inv(IV'*IV);

tmp1=x1'*IV;
t=(tmp1*miolo*x1)\(tmp1*miolo*y);
mvalold = exp(x1*t);
%oldt2 = zeros(size(theta2));
%mvalold = exp(mvalold);
save([data_dir 'mvalold.mat'],'mvalold');
%save([data_dir 'ps2.mat'],'s_jt','x1','x2','v');

clear mid y outshr t oldt2 mvalold temp sum1

%argums.invA=invA;
argums.theti=theti;
argums.thetj=thetj;
argums.x1=x1;
argums.x2=x2;
argums.IV=IV;
argums.ns=ns;
argums.vfull=vfull;
argums.cdindex=cdindex;
argums.cdid=cdid;
argums.data_dir=data_dir;
argums.s_jt=s_jt;


%tic
%[delta,its,norms]=tool_contr_PAR2(theta2,argums);
delta = meanval_PAR4(theta2,argums);
%toc
temp1 = x1'*IV;
theta1 = (temp1*miolo*x1)\(temp1*miolo*delta);
clear temp1 
gmmresid = delta - x1*theta1;
save([data_dir 'gmmresid.mat'],'gmmresid')

options = optimset('GradObj','on','DerivativeCheck','off','FinDiffType','central','TolFun',1e-6,'TolX',1e-6,'MaxFunEvals',1000);

%gmmobj(theta2,argums)

tic

[theta2,funcval2,exflag2,output2] = ktrlink(@(theta2) gmmobj(theta2,argums),theta2,[],[],[],[],[],[],[],options,[data_dir 'ktropts-a.txt']);

comp_t = toc/60;
disp(['Tempo de Execução: ' num2str(comp_t) ' minutos']);
% Carregando 
delta = meanval_PAR4(theta2,argums);
temp1 = x1'*IV;
temp2 = delta'*IV;


theta0 = (temp1*miolo*x1)\temp1*miolo*delta;

% Robust SE
gmmresid=delta-x1*theta0;

tic
miolo=inv(IV'*diag(gmmresid.^2)*IV);
disp(['Tempo de execução - direto: ' num2str(toc/60)]);
% 
% tic
% miolo2=wmatrix(IV,gmmresid,[1:size(IV,1)]');
% disp(['Tempo de execução - função: ' num2str(toc/60)]);
% 
% tic
% miolo3=wmatrix_mex(IV,gmmresid,[1:size(IV,1)]');
% disp(['Tempo de execução - função compilada: ' num2str(toc/60)]);


theta1=(temp1*miolo*temp1')\temp1*miolo*temp2';

argums.invA=miolo;
argums.v=v;
argums.perc=1;
argums.amount=0.01;
argums.demogr=[];
argums.dfull=dfull;
argums.theta1=theta1;

[aaa,tmp]=subst_blp_PAR(theta2,argums);
elasts2(:,mmm)=tmp;
elasts_pc=quantile(elasts2(:,1),[.10 .25 .50 .75 .90]);

vcov = var_cov(theta2,argums);
se = sqrt(diag(vcov));

temp=zeros(2*(size(theta1,1)+length(theta2)),1);
temp(1:2:end-1,:)=[theta1;theta2];
temp(2:2:end,:)=se;

coefsses=temp;

%scatter(elasts_IV,elasts2);title('Elasts LOGIT X Elasts RCL'); xlabel('LOGIT'); ylabel('RCL');
%scatter(x2(:,1),elasts_IV);title('Preço X Elasts LOGIT'); xlabel('Preço 1000SKR'); ylabel('Elasts Logit');
scatter(x2(:,1),elasts2);title('Preço X Elasts RCL'); xlabel('Preço 1000SKR'); ylabel('Elasts RCL');

