function [margeff,elast]=subst_blp2(theta2,args)
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
v=args.v;
demogr=args.demogr;



theta2w = full(sparse(theti,thetj,theta2));

mval=meanval_PAR4(theta2,args);
p = x2(:,1); 
temp1 = x1'*IV;
temp2 = mval'*IV;
theta1 = inv(temp1*invA*temp1')*temp1*invA*temp2';
expall = exp(mval*ones(1,ns));
% compute mu
mu=mufunc(theta2,args);
expmu=exp(mu);
eg=expall.*expmu;


temp = cumsum(eg);
sum1 = temp(cdindex,:);
sum1(2:size(sum1,1),:) = diff(sum1);
denom1 = 1./(1+sum1);
denom = denom1(cdid,:);
shari = eg.*denom;
sharsq=shari.^2;

margeff = theta1(end)+theta2(1).*vfull(:,1:ns);


%o1 = diag(ones(1,40)*margeff')*(shari*shari');

o1=sparse(size(shari,1),size(shari,1));
o2=sparse(size(shari,1),size(shari,1));

for i=min(cdid):max(cdid)
    ind=find(cdid==i);
    o1(ind,ind)=diag(sum(margeff(ind,:).*sharsq(ind,:),2));
    o2(ind,ind)=diag(sum(margeff(ind,:).*shari(ind,:),2));
end

%o1 = diag(ones(1,40)*margeff').*(shari*shari');
%o2 = diag(sum((margeff.*shari)'));

%o1 = spdiags((ones(1,ns)*margeff').*(shari*shari'),0,size(shari*shari'));
%o2 = spdiags(sum((margeff.*shari)'),0,size(shari*shari'));
omega = (o2-o1)/ns;

% Calculando as elasticidades
own_margeff=spdiags(omega,0);
elast=(own_margeff.*p)./s_jt;

