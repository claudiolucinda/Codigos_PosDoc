%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulations
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

load([data_dir 'dem_results.mat']);

simuls=dlmread([data_dir 'CFs_T_M.txt'],'\t',1,0);


Data.mg_cost=mg_cst;
Data.gmmresid=gmmresid;
Data.income=other2(:,1);
Data.nestid=nestid;
Data.pos=17;
Data.alpha2=alpha;
Data.ns=ns;

[p_init,s_init,checker]=Solver_BERT(Data.price,th12RCNL2S,Data.x1,Data.x2,Data);

[lsum_init]=Logsum_RCNL(p_init,th12RCNL2S,Data.x1,Data.x2,alpha,Data);

save([data_dir 'sim_results.mat']);

Data.tax=simuls(:,1);
[p_scen01,s_scen01,checker01]=Solver_BERT(Data.price,th12RCNL2S,Data.x1,Data.x2,Data);

[lsum_s_scen01]=Logsum_RCNL(p_init,th12RCNL2S,Data.x1,Data.x2,alpha,Data);

save([data_dir 'sim_results.mat']);

Data.tax=dlmread([data_dir 'tax_rate_T_M.txt'],'\t',1,0);
x2_alt=Data.x2;
x2_alt(:,1)=simuls(:,2);
[p_scen02,s_scen02,checker02]=Solver_BERT(Data.price,th12RCNL2S,Data.x1,x2_alt,Data);

[lsum_s_scen02]=Logsum_RCNL(p_init,th12RCNL2S,Data.x1,x2_alt,alpha,Data);

save([data_dir 'sim_results.mat']);

