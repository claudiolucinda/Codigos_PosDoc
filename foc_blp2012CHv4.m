function out=foc_blp2012CHv4(pr,theta2,args)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Código Base para o cálculo das matrizes de efeitos marginais
% Deve ser rodado depois da Estimação do BLP
% Copyright 2008 Cláudio R. Lucinda
% FGV/EESP e EAESP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs:
%   out - vetor de valores das CPO's por produto
% Inputs:
% theta2 - Vetor Final de Coeficientes
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
%   s_jt - vetor de participações de mercado
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


% Gerando o preço efetivamente enfrentado pelo consumidor
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
    




