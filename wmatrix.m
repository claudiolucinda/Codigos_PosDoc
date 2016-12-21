function w=wmatrix(IV,gmmresid,clustvar) %#codegen
%#codegen


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Função para cálculo mais rápido da matriz de ponderação
% A idéia é transformar isso em uma MEX - para rodar
% Inputs
% - IV - Matriz de Variáveis instrumentais
% - gmmresid - Vetor de Resíduos do GMM
% - cluster - vetor indicando os clusters - se quiser a versão robusta,
% colocar um vetor crescente por uma unidade com dim(linha)=dim(linha) de IV
% Outputs
% - w - Matriz de Ponderação
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coder.extrinsic('sparse','dummyvar')
temp1=max(max(clustvar,[],2),[],1);
temp2=size(IV,1);
if temp1==temp2
    w=inv(IV'*diag(gmmresid(:,1).^2)*IV);
    
else
    clust_dum=sparse(dummyvar(clustvar));
    pt_1=sparse((gmmresid*ones(1,size(IV,1)))*clust_dum');
    w=inv(IV'*(pt_1*pt_1')*IV);
end


