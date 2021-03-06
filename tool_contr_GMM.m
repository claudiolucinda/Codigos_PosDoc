    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arquivo Input-Output
% Copyright 2009 Cl�udio R. Lucinda
% FGV-EESP e EAESP
% Vers�o 13/03
% Agora puxando diretamente do Pharmaeco os dados
% Cen�rio Brasil vs. GCR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(RandStream('mcg16807', 'Seed',0));

clear; clc;
rand('seed',0);

cd 'C:\Users\claudiolucinda\Documents\Disc Taxes\trunk\'
% Se usando processamento paralelo descomentar aqui
matlabpool close force local
matlabpool open 2


% N�o chamando o gerador de draws
draw_dem=0;
if draw_dem==1
    draw_generator_jun2012
end

cd 'C:\Users\claudiolucinda\Documents\Disc Taxes\trunk\'


pharma_dir='Z:\Endog Chars\';
data_dir='C:\Users\claudiolucinda\Documents\Disc Taxes\trunk\';
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
ns=100;
nr_esp=1;
specs_wrk=[3 10];

y=data(:,1);
regrs=data(:,3:end);
x=[ones(size(y)) regrs];
x2=x(:,end); %x(:,2)  ones(size(x,1),1)];
x1=x(:,1:end);
s_jt=data2(:,1);


nobs=size(x,1);
ncoefs=size(x,2)+size(x2,2);

[n,K] = size(x2);

cdid=data2(:,2);
cdindex=find(diff(cdid));
cdindex=[cdindex;size(cdid,1)];



stored_draws=0;
if stored_draws==1
    v=load([pharma_dir2 'sigdraws2.mat']);
    
else
	vtmp=random('Normal',0,1,max(data2(:,2)),ns/2);
    v=[vtmp -vtmp];
    save([pharma_dir2 'sigdraws2.mat'],'v');
    clear vtmp
    
    
end
vfull=v(cdid,:);



%draws=load([data_dir 'draws_light4.mat']);
%d=[draws.d_draws_03;draws.d_draws_04;draws.d_draws_05;draws.d_draws_06;draws.d_draws_07;draws.d_draws_08;draws.d_draws_09];
dfull=[];
argums.dfull=dfull;


coefsses=zeros(2*ncoefs,nr_esp);
elasts2=zeros(nobs,nr_esp);

mmm=1;
kkk=1;
j=specs_wrk(kkk,1);
k=specs_wrk(kkk,2);
%z=[instvars(:,j) instvars(:,k) instvars(:,j).^2 instvars(:,k).^2 instvars(:,j).*instvars(:,k) instvars(:,end)];% instvars(:,end).*instvars(:,j)]; % instvars(:,end-1).*instvars(:,k)];
z=[instvars(:,j) instvars(:,k) instvars(:,j).*instvars(:,k) instvars(:,end)]; % instvars(:,j).*instvars(:,k) instvars(:,end)];% instvars(:,end).*instvars(:,j)]; % instvars(:,end-1).*instvars(:,k)];
% Truque para pegar e selecionar apenas as vari�veis que s�o LD

display(['Instrumentos S�o ' num2str(j) ' e ' num2str(k)]);

IV=[x1(:,1:end-1) z];
clear aaa


% Starting Values vari�vel RC - pre�o

theta2w=1e-5; % eps eps]';

[theti, thetj, theta2]=find(theta2w);

% Starting Values - Vari�veis n�o RC
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
elasts_IV=theta1_IV(end)*(1-s_jt).*x2;
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
teste=1;
results=zeros(100,3);
for thetatemp=1e-10:1e-4:1e-2
    results(teste,1)=thetatemp;
    [ftemp, df]=gmmobj(thetatemp,argums);
    results(teste,2)=ftemp;
    results(teste,3)=df;
    teste=teste+1;
    
end

scatter(results(:,1),results(:,2)); hold on
scatter(results(:,1),results(:,3)); hold off
