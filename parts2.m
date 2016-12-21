function share=parts2(theta2,args)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Código Base para o cálculo das matrizes de efeitos marginais
% Deve ser rodado depois da Estimação do BLP
% Copyright 2008 Cláudio R. Lucinda
% FGV/EESP e EAESP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs:
%   margeff - matriz de efeitos marginais do preço na utilidade indireta
%   condicional (alphas)
%   elast - vetor de elasticidades-preço próprias
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