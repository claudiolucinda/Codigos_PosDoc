function w=wmatrix(IV,gmmresid,clustvar) %#codegen
%#codegen


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fun��o para c�lculo mais r�pido da matriz de pondera��o
% A id�ia � transformar isso em uma MEX - para rodar
% Inputs
% - IV - Matriz de Vari�veis instrumentais
% - gmmresid - Vetor de Res�duos do GMM
% - cluster - vetor indicando os clusters - se quiser a vers�o robusta,
% colocar um vetor crescente por uma unidade com dim(linha)=dim(linha) de IV
% Outputs
% - w - Matriz de Pondera��o
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


