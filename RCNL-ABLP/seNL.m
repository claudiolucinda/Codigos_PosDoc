function se     = seNL(parameters,Data)

warning('off','MATLAB:nearlySingularMatrix');

%global args

fnames=fieldnames(Data);
for i=1:length(fnames)
    eval([fnames{i} '=Data.' fnames{i} ';']);
end



load([data_dir 'mvaloldNL.mat'],'meanval_FINAL');
load([data_dir 'gmmresid.mat'],'gmmresid');

nx1=size(x1,2);
nx2=size(x2,2);
theta2=parameters(nx1+1:nx1+nx2-1,:);

Z = size(IV,2);
temp = jacobNL(meanval_FINAL,theta2,Data);
a = [x1 temp]'*IV;

IVres = IV.*(gmmresid*ones(1,Z));
b = IVres'*IVres;

% Original
varcovar = (a*W*a')\a*W*b*W*a'/(a*W*a');
se       = sqrt(diag(varcovar));

end