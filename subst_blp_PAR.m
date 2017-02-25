function [margeff,elast]=subst_blp_PAR(theta2,args)
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
warning('off','MATLAB:nearlySingularMatrix');

%global args

fnames=fieldnames(args);
for i=1:length(fnames)
    eval([fnames{i} '=args.' fnames{i} ';']);
end



theta2w = full(sparse(theti,thetj,theta2));

mval=meanval_PAR4(theta2,args);
p = x1(:,end); 
miolo=miolo2*IV';
temp1 = x1'*IV;
temp2 = mval'*IV;
theta1 = (temp1*miolo*x1)\temp1*miolo*mval;
expall = exp(mval*ones(1,ns));
% compute mu
mu=mufunc(theta2,args);
expmu=exp(mu);
eg=expall.*expmu;
temp_oo=dummyvar(cdid);
sumexp=temp_oo'*eg;
denom = 1./(1+sumexp);
sum1=temp_oo*denom;
shari = eg.*sum1;
%sharsq=shari.^2;
margeff=theta1(end);



%o1 = diag(ones(1,40)*margeff')*(shari*shari');

o1=sparse(size(shari,1),size(shari,1));
o2=sparse(size(shari,1),size(shari,1));

for i=min(cdid):max(cdid)
   mktsubs{i}=cdid==i;
   meff{i}=margeff;
   sqs{i}=shari(mktsubs{i},:);
    
end

parfor i=min(cdid):max(cdid)
    meffs=meff{i};
    shars=sqs{i};
    o1_tmp=(meffs.*shars)*shars';
    o2_tmp=diag(sum(meffs.*shars,2));
    omg{i}=num2cell(o2_tmp-o1_tmp);
end

omega=sparse(size(shari,1),size(shari,1));

for i=min(cdid):max(cdid)
    omega(mktsubs{i},mktsubs{i})=cell2mat(omg{i})/ns;
end


%o1 = diag(ones(1,40)*margeff').*(shari*shari');
%o2 = diag(sum((margeff.*shari)'));

%o1 = spdiags((ones(1,ns)*margeff').*(shari*shari'),0,size(shari*shari'));
%o2 = spdiags(sum((margeff.*shari)'),0,size(shari*shari'));
%omega = (o2-o1)/ns;

% Calculando as elasticidades
own_margeff=spdiags(omega,0);
elast=(own_margeff.*p)./s_jt;

