    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arquivo Input-Output
% Copyright 2009 Cláudio R. Lucinda
% FGV-EESP e EAESP
% Versão 13/03
% Agora puxando diretamente do Pharmaeco os dados
% Cenário Brasil vs. GCR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(RandStream('mcg16807', 'Seed',0));
% Se usando processamento paralelo descomentar aqui
matlabpool close force local
matlabpool open 2

clear; clc;
rand('seed',0);

cd  'C:\Users\claudiolucinda\Documents\Disc Taxes\branches\branch06\'
data_dir='C:\Users\claudiolucinda\Documents\Disc Taxes\branches\branch06\';



% Não chamando o gerador de draws
draw_dem=1;
if draw_dem==1
    demogr=load([data_dir 'draws_dem.mat']);
    demogr=demogr.draws_fin;
else 
    'draw_gen_2013.m';
    clear; clc;
    demogr=load([data_dir 'draws_dem.mat']);
    demogr=demogr.draws_fin;
end




pharma_dir='Z:\Endog Chars\';
path(data_dir,path);


pharma_dir2='Z:\Endog Chars\BLP_Intcodes\';

data=dlmread([pharma_dir 'modelvars.txt'],'\t',1,0);
marca=dlmread([pharma_dir 'brand_dummies.txt'],'\t',1,0);
data2=dlmread([pharma_dir 'otherdata.txt'],'\t',1,0);
instvars=dlmread([pharma_dir 'instvars.txt'],'\t',1,0);

% Importando as coisas
% data_alt=dlmread([pharma_dir 'outdataSWE2012alt.txt'],'\t',1,0);
% instvars_alt=dlmread([pharma_dir 'instmat2012alt.txt'],'\t',1,0);
% data_alt2=dlmread([data_dir 'outdata2SWE2012alt.txt'],'\t',1,0);
% tfuel_vars=dlmread([pharma_dir 'typefuel.txt'],'\t',1,0);
% summer=dlmread([pharma_dir 'summer.txt'],'\t',1,0);
% ksks=dlmread([pharma_dir 'ks_allb.txt'],'\t',1,0);
% otherstuff=dlmread([pharma_dir 'otherstuffb.txt'],'\t',1,0);
% marca=dlmread([pharma_dir 'marca.txt'],'\t',1,0);
% contrafact=dlmread([pharma_dir 'contrafact.txt'],'\t',1,0);
% data=dlmread([pharma_dir 'outdataSWE2012h.txt'],'\t',1,0);
% instvars=dlmread([pharma_dir 'instmat2012h.txt'],'\t',1,0);
% data2=dlmread([pharma_dir 'outdata2SWE2012h.txt'],'\t',1,0);
% trend=dlmread([pharma_dir 'timetrend.txt'],'\t',1,0);
% contrafact2=dlmread([pharma_dir 'contrafact2.txt'],'\t',1,0);


dum_marca=dummyvar(marca);
ns=250;
nr_esp=1;
specs_wrk=[3 10];

y=data(:,1);
regrs=data(:,3:end);
x=[ones(size(y)) regrs];
x2=[x(:,end) x(:,2)]; %x(:,2)  ones(size(x,1),1)];
x1=x(:,1:end-1);
s_jt=data2(:,1);

%x2(:,1)=x2(:,1)/1000;
%x1(:,end)=x1(:,end)/1000;

nobs=size(x,1);
ncoefs=size(x,2)+size(x2,2);

[n,K] = size(x2);

cdid=data2(:,2);
cdindex=find(diff(cdid));
cdindex=[cdindex;size(cdid,1)];



stored_draws=1;
if stored_draws==1
    v=load([pharma_dir2 'sigdraws2.mat']);
    v=v.v;
else
	vtmp=random('Normal',0,1,max(data2(:,2)),ns/2);
    v=[vtmp -vtmp];
    save([pharma_dir2 'sigdraws2.mat'],'v');
    clear vtmp
end

vfull=v(cdid,:);
dfull=demogr(cdid,:);
argums.dfull=dfull;


coefsses=zeros(2*ncoefs,nr_esp);
elasts2=zeros(nobs,nr_esp);

mmm=1;
kkk=1;
j=specs_wrk(kkk,1);
k=specs_wrk(kkk,2);
z=[instvars(:,j) instvars(:,k) instvars(:,j).*instvars(:,k) instvars(:,end)]; % instvars(:,j).*instvars(:,k) instvars(:,end)];% instvars(:,end).*instvars(:,j)]; % instvars(:,end-1).*instvars(:,k)];

display(['Instrumentos São ' num2str(j) ' e ' num2str(k)]);

IV=[x1(:,1:end-1) z];
clear aaa


% Starting Values variável RC - preço

theta2w=[1 1e-5]'; % eps eps]';

[theti, thetj, theta2]=find(theta2w);

% Starting Values - Variáveis não RC
miolo=pinv(IV);

%invA = inv(IV'*IV);

tmp1=x1'*IV;
t=(tmp1*miolo*x1)\(tmp1*miolo*y);
mvalold = exp(x1*t);
%oldt2 = zeros(size(theta2));
%mvalold = exp(mvalold);
save([pharma_dir2 'mvalold.mat'],'mvalold');
%save([data_dir 'ps2.mat'],'s_jt','x1','x2','v');
theta1_IV=t;
elasts_IV=theta1_IV(end)*(1-s_jt).*(x2(:,1));
elasts_pc_IV=quantile(elasts_IV,[.1 .25 .50 .75 .90]);

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
argums.data_dir=pharma_dir2;
argums.s_jt=s_jt;


%tic
%[delta,its,norms]=tool_contr_PAR2(theta2,argums);
delta = meanval_PAR4(theta2,argums);
%toc
%delta = meanval(theta2,argums);
temp1 = x1'*IV;
%temp2 = delta'*IV;
theta1 = (temp1*miolo*x1)\(temp1*miolo*delta);
clear temp1 
gmmresid = delta - x1*theta1;
save([pharma_dir2 'gmmresid.mat'],'gmmresid')

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

