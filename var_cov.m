function f = var_cov(theta2,args)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Esta Matriz Calcula a Matriz Vari�ncia-Covari�ncia das EStimativas
% % Adaptada e otimizada por Cl�udio R. Lucinda
% FGV-EESP e EAESP
% 2008/09
% Um Output:
% f - Matriz VC
% df - Gradiente da Fun��o - Para otimizadores que precisam
%
% Inputs:
% theta2 - Vetor de Parametros (o theta2w "amassado")
% args - estrutura com os seguinte campos
%   invA - Matriz de Pondera��o do GMM (Weight Matrix)
%   theti - ind�ce linha dos termos iguais a zero no theta2w
%   thetj - �ndice coluna dos termos iguais a zero no theta2w
%   x1 - Matriz de Dados da parte linear da estima��o
%   x2 - Matriz de Dados da parte n�o-linear da estima��o
%   IV - instrumentos
%   ns - n�mero de indiv�duos simulados
%   vfull - parte n�o observ�vel j� expandida
%   dfull - parte draws demografia j� expandida
%   cdindex - vetor que quebra os mercados
%   cdid - vetor que diz a que mercado pertence cada observa��o
%   data_dir - string com o diret�rio aonde est�o os dados
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

