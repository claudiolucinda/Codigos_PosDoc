function f = jacob(mval,theta2,args)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Esta função calcula o valor do Jacobiano da Função GMM
% Adaptada e otimizada por Cláudio R. Lucinda
% FGV-EESP e EAESP
% 2008/09
% Um Output:
% f - matriz Jacobiana
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

%invA=args.invA;
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

[n,K]=size(x2);
theta2w = full(sparse(theti,thetj,theta2));
mu=mufunc(theta2,args);
expmu = exp(mu);

eg = expmu.*(exp(mval)*ones(1,ns));
temp_oo=dummyvar(cdid);
sumexp=temp_oo'*eg;
denom = 1./(1+sumexp);
sum1=temp_oo*denom;
shares = eg.*sum1;
clear expmu

J = size(theta2w,2) - 1;
f1 = zeros(size(cdid,1),K*(J + 1));

% computing (partial share)/(partial sigma)


if K>1
for i = 2:K
    l=i-1;
    xv = (x2(:,i)*ones(1,ns)).*vfull(:,(l-1)*ns+1:l*ns);    
    sumxv=temp_oo'*(xv.*shares);
    sumxvexp=temp_oo*sumxv;
    f1(:,i)=mean((shares.*(xv-sumxvexp)),2)';
    clear xv temp sum1
end
end

j=1;
yp=dfull(:,(j-1)*ns+1:j*ns)-(x2(:,j)*ones(1,ns));
neg=yp<=0;
yp=yp.*(1-neg)+neg;
yp=(log(yp)-log(dfull(:,(j-1)*ns+1:j*ns)));
dtmp=yp.*(1-neg);
sumxd=temp_oo'*(dtmp.*shares);
sumxdexp=temp_oo*sumxd;
% dtmp=dfull(:,(j-1)*ns+1:j*ns);
% xd=(x2(:,j)*ones(1,ns)).*dtmp;
% sumxd=temp_oo'*(xd.*shares);
% sumxdexp=temp_oo*sumxd;
f1(:,j)=mean((shares.*(dtmp-sumxdexp)),2)';

% temp=cumsum(xd.*shares);
% sum1 = temp(cdindex,:);
% sum1(2:size(sum1,1),:) = diff(sum1);
% f1(:,j) = mean((shares.*(xd-sum1(cdid,:))),2)';
% 
% 
% i=1;
% xv = (x2(:,i)*ones(1,ns)).*vfull(:,1:ns);
% sumxv=temp_oo'*(xv.*shares);
% sumxvexp=temp_oo*sumxv;
% f1(:,i)=mean((shares.*(xv-sumxvexp)),2)';
%temp = cumsum(xv.*shares);
%sum1 = temp(cdindex,:);
%sum1(2:size(sum1,1),:) = diff(sum1);
%f1(:,i) = mean((shares.*(xv-sum1(cdid,:))),2)';
clear xd temp sum1



% If no demogr comment out the next para
% computing (partial share)/(partial pi)
% for j = 1:J
% d = dfull(:,ns*(j-1)+1:ns*j);    
% 	temp1 = zeros(size(cdid,1),K);
% 	for i = 1:K
% 		xd=(x2(:,i)*ones(1,ns)).*d;    
% 		temp = cumsum(xd.*shares);
% 		sum1 = temp(cdindex,:);
% 		sum1(2:size(sum1,1),:) = diff(sum1);
% 		temp1(:,i) = mean((shares.*(xd-sum1(cdid,:)))')';
% 		clear xd temp sum1
% 	end
% 	f1(:,K*j+1:K*(j+1)) = temp1;
% 	clear temp1
% end

rel = theti + (thetj - 1) * max(theti) ;

% computing (partial delta)/(partial theta2)

f = zeros(size(cdid,1),size(rel,1));
n = 1;
for i = 1:size(cdindex,1)
	temp = shares(n:cdindex(i),:);
	H1 = temp*temp';
	H = (diag(sum(temp,2)) - H1)/ns;
	f(n:cdindex(i),:) = - H\f1(n:cdindex(i),rel);
	n = cdindex(i) + 1;
end
