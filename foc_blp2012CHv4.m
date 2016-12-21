function out=foc_blp2012CHv4(pr,theta2,args)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C�digo Base para o c�lculo das matrizes de efeitos marginais
% Deve ser rodado depois da Estima��o do BLP
% Copyright 2008 Cl�udio R. Lucinda
% FGV/EESP e EAESP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs:
%   out - vetor de valores das CPO's por produto
% Inputs:
% theta2 - Vetor Final de Coeficientes
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
%   s_jt - vetor de participa��es de mercado
%   v - Vetor de random draws
%   demogr - Vetor de draws de demografia
%   dum_extra - matriz de brand dummies
%   CMg - Vetor de Custos Marginais
%   rebate - Vetor de Rebates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Desempacotando argums
fnames=fieldnames(args);
for i=1:length(fnames)
    eval([fnames{i} '=args.' fnames{i} ';']);
end


% Gerando o pre�o efetivamente enfrentado pelo consumidor
pcons=pr+rebate;

x1(:,end)=pcons;

%theta2w = full(sparse(theti,thetj,theta2));

mval1=x1*theta1;
expall = exp(mval1*ones(1,ns));

margeff = theta1(end)+theta2(1)*vfull;

args.x2(:,1)=pcons;
mu=mufunc(theta2,args);

eg = expall.*exp(mu);
sum1 = sum(eg,1);
denom = 1./(1+ones(size(eg,1),1)*sum1);
shari = eg.*denom;


o2=diag(diag(margeff*shari'));
o1=(margeff.*shari)*shari';
omega=(o2-o1)/ns;
own_mat1=dum_extra*dum_extra';
mgef=omega.*own_mat1;
out=mgef*(pr-cmg)+mean(shari,2);
    




