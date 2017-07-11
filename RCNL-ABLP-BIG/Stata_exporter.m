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

load([data_dir 'sim_results.mat']);

key=dlmread([data_dir 'key.txt'],'\t',1,0);

data_exp=[key p_init s_init p_scen01 s_scen01 p_scen02 s_scen02];

lsum_exp=[lsum_init lsum_s_scen01 lsum_s_scen02];

dlmwrite([data_dir 'sim_results.txt'],data_exp,'\t');

dlmwrite([data_dir 'logsums.txt'],lsum_exp,'\t');
