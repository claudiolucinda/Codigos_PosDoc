function [mkup,cmg]=supply2012CH(theta2,args)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C�digo Base para o c�lculo das matrizes de efeitos marginais
% Deve ser rodado depois da Estima��o do BLP
% Copyright 2008 Cl�udio R. Lucinda
% FGV/EESP e EAESP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs:
%   mkup - Vetor de Markups
%   cmg - Vetor de Custos Marginais
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
%   dum_extra - matriz de brand dummies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


invA=args.invA;
theti=args.theti;
thetj=args.thetj;
x1=args.x1;
x2=args.x2;
IV=args.IV;
ns=args.ns;
vfull=args.vfull;
dfull=args.dfull;
cdindex=args.cdindex;
cdid=args.cdid;
data_dir=args.data_dir;
s_jt=args.s_jt;
%v=args.v;
%demogr=args.demogr;
dum_extra=args.dum_extra;
theta1=args.theta1;



%theta2w = full(sparse(theti,thetj,theta2));

%mval=meanval(theta2,args);
p=x2(:,1);
mval1=x1*theta1;
expall = exp(mval1*ones(1,ns));

margeff = theta1(end)+theta2(1)*vfull;

mu=mufunc(theta2,args);

eg = expall.*exp(mu);

temp = cumsum(eg);
sum1 = temp(cdindex,:);
sum1(2:size(sum1,1),:) = diff(sum1);
denom1 = 1./(1+sum1);
denom = denom1(cdid,:);
shari = eg.*denom;
%sharsq=shari.^2;


% Fazendo a matriz aqui ano por ano

mkup_BRL=zeros(size(shari,1),1);
cmg=zeros(size(mkup_BRL));
mkup=zeros(size(mkup_BRL));

for i=min(cdid):max(cdid)
    ind=find(cdid==i);
    o2=diag(diag(margeff(ind,:)*shari(ind,:)'));
    o1=(margeff(ind,:).*shari(ind,:))*shari(ind,:)';
    omega=(o2-o1)/ns;
    own_mat1=dum_extra(ind,:)*dum_extra(ind,:)';
    mkup_BRL(ind,:)=-inv(own_mat1.*omega)*(mean(shari(ind,:),2));
    mkup(ind,:)=mkup_BRL(ind,:)./p(ind,:);
    cmg(ind,:)=p(ind,:)-mkup_BRL(ind,:);
end
