function f = deltaNL_B(thetaNL,Data)

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
k   = 100;
km  = 1e-14;
i   = 0;
if max(isnan(delta0NL)|isinf(delta0NL))  == 1
    deltastart           =  zeros(Data.nobs,1);
else
    deltastart           = delta0NL;
end

% Unpack
logobsshare              = log(sj);

while k > km
    % Market Share
    [~,~,sh,~,~] = NLShareCalculation(thetaNL,deltastart,Data);
    %sh      = reshape((sum(sh,2)),[Data.nobs 1]);   % 2D market shares
    % Modified contraction mapping for the random coefficient nested logit
    % theta(1) is sigseg
    delta1  = deltastart+((1-max(thetaNL(end)))*(logobsshare-log(sh)));
    if max(isnan(delta1)) == 1
        disp('No Convergence - delta calculation failed: OVERFLOW')
        break
    end
    i = i + 1;
    if i>2500
        disp('No Convergence - delta convergence failed')
        break
    end
    k           = max(abs(delta1-deltastart));
    deltastart  = delta1;
end

disp(['# of iterations for delta convergence:  ' num2str(i)])

delta0NL = delta1;
save([data_dir 'mvaloldNL.mat'],'delta0NL');

f=delta1;