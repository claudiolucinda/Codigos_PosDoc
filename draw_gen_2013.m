    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arquivo Input-Output
% Copyright 2013 Cláudio R. Lucinda
% FEA-RP/USP
% Versão 01 - Apr. 2013
% Agora puxando diretamente do Pharmaeco os dados
% Cenário Brasil vs. GCR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(RandStream('mcg16807', 'Seed',0));

clear; clc;
rand('seed',0);

cd  'C:\Users\claudiolucinda\Documents\Disc Taxes\branches\branch06\'
%cd 'C:\Users\claudiolucinda\Documents\Disc Taxes\trunk\'
% Se usando processamento paralelo descomentar aqui
matlabpool close force local
matlabpool open 2


data_dir='C:\Users\claudiolucinda\Documents\Disc Taxes\branches\branch06\';
path(data_dir,path);

ns=250;
params=dlmread([data_dir 'draws_guideline.txt'],'\t',1,0);

lims_inf=[1 20 40 60 80 100 120 140 160 180 200 220 240 260 280 300 320 340 360 380 400 500 600 800 1000]';
lims_sup=[lims_inf(2:end)-1;6000];

data=params(:,3:27);

draws_fin=zeros(size(data,1),2*ns);
%i=3;


parfor i=1:size(data,1)
temprow=data(i,:);
dtemp=[];
for j=1:size(data,2);
    % Versão com Antithetics
    dunif11=random('unif',0,1,1,data(i,j));
    dunif12=1-dunif11;
    dtemp11=lims_inf(j,1)+dunif11*(lims_sup(j,1)-lims_inf(j,1)-eps);
    dtemp12=lims_inf(j,1)+dunif12*(lims_sup(j,1)-lims_inf(j,1)-eps);
    %dtemp2=random('unif',lims_inf(j,1),lims_sup(j,1)-eps,1,data(i,j));
    dtemp=[dtemp dtemp11 dtemp12];
end
if size(dtemp,2)>2*ns
    sel = randperm(2*ns);
    sel = sel(1:size(dtemp,2)-2*ns);
    dtemp(:,sel)=[];
end
draws_fin(i,:)=dtemp;
   
end
sel=randperm(2*ns);
sel=sel(1:ns);
draws_fin(:,sel)=[];
draws_fin=draws_fin;
save([data_dir 'draws_dem.mat'],'draws_fin');