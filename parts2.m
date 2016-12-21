function share=parts2(theta2,args)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C�digo Base para o c�lculo das matrizes de efeitos marginais
% Deve ser rodado depois da Estima��o do BLP
% Copyright 2008 Cl�udio R. Lucinda
% FGV/EESP e EAESP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs:
%   margeff - matriz de efeitos marginais do pre�o na utilidade indireta
%   condicional (alphas)
%   elast - vetor de elasticidades-pre�o pr�prias
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%invA=args.invA;
theti=args.theti;
thetj=args.thetj;
x1=args.x1;
x2=args.x2;
%IV=args.IV;
ns=args.ns;
vfull=args.vfull;
dfull=args.dfull;
cdindex=args.cdindex;
cdid=args.cdid;
data_dir=args.data_dir;
%s_jt=args.s_jt;
v=args.v;
demogr=args.demogr;
theta1=args.theta1;



%theta2w = theta2;
theta2w = full(sparse(theti,thetj,theta2));

mval=x1*theta1;
mu=mufunc(theta2,args);
%[n k] = size(x2);
%j = size(theta2w,2)-1;
mval1 = mval*ones(1,ns) + mu;
eg = exp(mval1);
% compute individual price sensitivity, alphai (will be used below)

%v_long = reshape(v',size(cdindex,1)*ns,1);
%margeff = theta2(1)./dfull;
%alphai = reshape(alphai,ns, size(cdindex,1))';

%eg = expall.*exp(p*ones(1,ns).*margeff);
temp = cumsum(eg);
sum1 = temp(cdindex,:);
sum1(2:size(sum1,1),:) = diff(sum1);
denom1 = 1./(1+sum1);
denom = denom1(cdid,:);
shari = eg.*denom;
share=sum(shari,2)/ns;