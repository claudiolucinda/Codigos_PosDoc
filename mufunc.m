function mu=mufunc(theta2,args)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Esta fun��o calcula a parte mu da utilidade m�dia
% Coloquei esta fun��o aqui para facilitar o processo
% Se tiver que mexer depois, � em um lugar s�
% Output:
%   mu=uma matriz com (size(x2,1)) linhas e (ns) colunas
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
%   vfull - parte n�o observ�vel j� expandida - ela � (size(x2,1)) linhas e
%   size(x2,2)*ns colunas (DIFEREN�A AQUI)
%   dfull - parte draws demografia j� expandida
%   cdindex - vetor que quebra os mercados
%   cdid - vetor que diz a que mercado pertence cada observa��o
%   data_dir - string com o diret�rio aonde est�o os dados
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fnames=fieldnames(args);
for i=1:length(fnames)
    eval([fnames{i} '=args.' fnames{i} ';']);
end

%[n,K] = size(x2);
%j = size(theta2w,2)-1;
[n,K] = size(x2);
%j = size(theta2w,2)-1;
if isempty(dfull)==0
    j=1;
    mu=((x2(:,j)*ones(1,ns)).*dfull(:,(j-1)*ns+1:j*ns)*theta2(j));

    for j=2:K
        mu=mu+((x2(:,j)*ones(1,ns)).*vfull(:,(j-2)*ns+1:(j-1)*ns)*theta2(j));

    end
else
    mu=0;
    for j=1:K
        mu=mu+((x2(:,j)*ones(1,ns)).*vfull(:,(j-1)*ns+1:j*ns)*theta2(j));

    end
end
    
end


