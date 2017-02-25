function f = var_cov(theta2,args)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Esta Matriz Calcula a Matriz Variância-Covariância das EStimativas
% % Adaptada e otimizada por Cláudio R. Lucinda
% FGV-EESP e EAESP
% 2008/09
% Um Output:
% f - Matriz VC
% df - Gradiente da Função - Para otimizadores que precisam
%
% Inputs:
% theta2 - Vetor de Parametros (o theta2w "amassado")
% args - estrutura com os seguinte campos
%   invA - Matriz de Ponderação do GMM (Weight Matrix)
%   theti - indíce linha dos termos iguais a zero no theta2w
%   thetj - índice coluna dos termos iguais a zero no theta2w
%   x1 - Matriz de Dados da parte linear da estimação
%   x2 - Matriz de Dados da parte não-linear da estimação
%   IV - instrumentos
%   ns - número de indivíduos simulados
%   vfull - parte não observável já expandida
%   dfull - parte draws demografia já expandida
%   cdindex - vetor que quebra os mercados
%   cdid - vetor que diz a que mercado pertence cada observação
%   data_dir - string com o diretório aonde estão os dados
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','MATLAB:nearlySingularMatrix');

%global args

fnames=fieldnames(args);
for i=1:length(fnames)
    eval([fnames{i} '=args.' fnames{i} ';']);
end



load([data_dir 'mvalold.mat'])
load([data_dir 'gmmresid.mat'])
%N = size(x1,1);
Z = size(IV,2);
temp = jacob(mvalold,theta2,args);
a = [x1 temp]'*IV;
IVres = IV.*(gmmresid*ones(1,Z));
b = IVres'*IVres;

% Original
f = (a*invA*a')\a*invA*b*invA*a'/(a*invA*a');

%f = (a*miolo*[x1 temp])\a*invA*b*invA*a'/inv(a*miolo*[x1 temp]);

