function g = jacob_solvopt(theta2)
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

global args

fnames=fieldnames(args);
for i=1:length(fnames)
    eval([fnames{i} '=args.' fnames{i} ';']);
end

%p=x2(:,1);

mval1 = load([data_dir '\mvalold.mat']);
mval=mval1.mvalold;

K=size(x2,2);
theta2w = full(sparse(theti,thetj,theta2));
mu=mufunc(theta2);
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
%f2 = zeros(size(cdid,1),K*(J + 1));

%if optimalinstrument==0;
    
if isempty(dfull)==0;
    for i = 2:K
        l=i-1;
        xv = (x2(:,i)*ones(1,ns)).*vfull(:,(l-1)*ns+1:l*ns);
        sumxv=temp_oo'*(xv.*shares);
        sumxvexp=temp_oo*sumxv;
        f1(:,i)=mean((shares.*(xv-sumxvexp)),2)';
        %clear xv temp sum1
        %f2(:,i)=mean((((p-cmg)./cmg)*ones(1,ns)).*((1./mu).*vfull(:,1:ns)+(1./(1-BigTheta*shares)).*(BigTheta*(shares.*(xv-sumxvexp)))),2);
    end
    clear xv temp sum1

    j=1;
    
    xd=(x2(:,j)*ones(1,ns)).*dfull(:,(j-1)*ns+1:j*ns);
    sumxd=temp_oo'*(xd.*shares);
    sumxdexp=temp_oo*sumxd;
    f1(:,j)=mean((shares.*(xd-sumxdexp)),2)';
else
    for i = 1:K
        xv = (x2(:,i)*ones(1,ns)).*vfull(:,(K-1)*ns+1:K*ns);
        sumxv=temp_oo'*(xv.*shares);
        sumxvexp=temp_oo*sumxv;
        f1(:,i)=mean((shares.*(xv-sumxvexp)),2)';
    end
    clear xv temp sum1
end    
    
    
    
    
    
    
    
    
    
% elseif optimalinstrument==1;
% 
%     if K>1
%         for i = 2:K
%             l=i-1;
%             xvtemp = xv(:,(l-1)*ns+1:l*ns);
%             sumxv=temp_oo'*(xvtemp.*shares);
%             sumxvexp=temp_oo*sumxv;
%             f1(:,i)=mean((shares.*(xvtemp-sumxvexp)),2)';
%             %clear xv temp sum1
%             %f2(:,i)=mean((((p-cmg)./cmg)*ones(1,ns)).*((1./mu).*vfull(:,1:ns)+(1./(1-BigTheta*shares)).*(BigTheta*(shares.*(xv-sumxvexp)))),2);
%         end
%     end
%     
%     clear xv temp sum1
%     
%     j=1;
%     sumxd=temp_oo'*(xd.*shares);
%     sumxdexp=temp_oo*sumxd;
%     f1(:,j)=mean((shares.*(xd-sumxdexp)),2)';
%     
%     
%     

%f2(:,j)=mean((((p-cmg)./cmg)*ones(1,ns)).*((1./mu).*dtmp+(1./(1-BigTheta*shares)).*(BigTheta*(shares.*(xd-sumxdexp)))),2);

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

g=f;
% Fazendo as coisas
% Acertar aqui as colunas pros diferentes coeficientes
% Miolo =(1./mu).*vfull+(1./(1-dummyset'*share)).*(shares.*(xv-sumxvexp))
% grad_supp=(((p-cmg)/cmg)*ones(1,ns)).*Miolo;
% whole=[f;grad_supp];
% 