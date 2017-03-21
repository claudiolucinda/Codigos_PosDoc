function f = deltaNL(thetaNL,Data)

% DeltaNL - Contraction Mapping for Random Coefficient Nested Logit
% Works like Logit apart for the weighting of the difference between
% calculated and observed shares
% Syntax:  f = deltaNL(thetaNL,Data)
%
% Inputs:
%    input1 - Non linear parameters
%    input2 - Data
%
% Outputs:
%    Mean Value Delta
%
% MAT-files required: mvaloldNL
% Subfunctions: NLShareCalculation
%

% Author: Laura Grigolon and Frank Verboven
% August 2012;


fnames=fieldnames(Data);
for i=1:length(fnames)
    eval([fnames{i} '=Data.' fnames{i} ';']);
end

load([data_dir 'mvaloldNL.mat'],'delta0NL');
if max(isnan(delta0NL)|isinf(delta0NL))  == 1
    deltastart           =  zeros(size(delta0NL));
else
    deltastart           = delta0NL;
end

tol=1e-13;

logobsshare              = log(sj);


warning off
i2=0;
norm=1;
while norm > tol && i2<2500
    i2 = i2+1; 
    [~,~,sh,~,~] = NLShareCalculation(thetaNL,deltastart,Data);
    delta1  = deltastart+((1-max(thetaNL(end)))*(logobsshare-log(sh)));
    
    norm = max(abs(delta1-deltastart));
    deltastart  = delta1;

end
disp(['Número de Iterações: ' num2str(i2)]);
warning on

delta0NL = delta1;
save([data_dir 'mvaloldNL.mat'],'delta0NL');

f=delta1;