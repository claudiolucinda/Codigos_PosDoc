function [se,varcovar]     = seNL(parameters,Data)

warning('off','MATLAB:nearlySingularMatrix');

%global args

fnames=fieldnames(Data);
for i=1:length(fnames)
    eval([fnames{i} '=Data.' fnames{i} ';']);
end



%load([data_dir 'mvaloldNL.mat'],'delta0NL');
%load([data_dir 'gmmresid.mat'],'gmmresid');

nx1=size(x1,2);
nx2=size(x2,2);
theta2=parameters(nx1+1:nx1+nx2,:);
deltaopt = deltaNL_PAR(theta2,Data);
gmmresid = deltaopt - x1*parameters(1:nx1,:);
%Z = size(IV,2);
temp = jacobNL(deltaopt,theta2,Data);
a = [-x1 temp]'*IV;

% residdiag = spdiags(gmmresid.*gmmresid ,0,length(gmmresid),length(gmmresid));
% sqrtresiddiagIV = (residdiag.^.5)*sparse(IV);
% [~,R] = qr(sqrtresiddiagIV,0);
% covg=R'*R;
% 
% covg     = zeros(size(Data.IV,2));
% for ii   = 1:length(Data.IV),
%      covg = covg + Data.IV(ii,:)'*Data.IV(ii,:)*(gmmresid(ii)^2);
% end

IVres = IV.*(gmmresid*ones(1,size(IV,2)));
covg = IVres'*IVres;

% Original
varcovar = (a*W*a')\a*W*covg*W*a'/(a*W*a');
se       = sqrt(diag(varcovar));

end