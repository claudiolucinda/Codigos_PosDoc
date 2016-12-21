function [f,df] = gmmobj(theta2,args)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Esta função calcula o valor da função objetivo GMM
% Adaptada e otimizada por Cláudio R. Lucinda
% FGV-EESP e EAESP
% 2008/09
% Dois Outputs:
% f - Valor da função objetivo
% df - Gradiente da Função - Para otimizadores que precisam
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
%   s_jt - vetor de participações de mercado
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

% Carregando 
delta = meanval_PAR4(theta2,args);

% O que fazer quando o algoritmo sai fora
if max(isnan(delta)) == 1
	f = 1e+10;
    gmmresid=(1e+10)*ones(size(delta));
    df=((1e+10)*ones(size(theta2)));         
    save([data_dir 'gmmresid.mat'],'gmmresid')
	   
else
    miolo=pinv(IV);
    temp1 = x1'*IV;
	theta1 = (temp1*miolo*x1)\temp1*miolo*delta;
    clear temp1 
	gmmresid = delta - x1*theta1;
	temp1 = gmmresid'*IV;
    f = temp1*miolo*gmmresid;
	%f = temp1*inv(IV'*IV)*temp1';
    clear temp1
	save([data_dir 'gmmresid.mat'],'gmmresid')
end

disp(['GMM:  ' num2str(f)])

% % Parte do gradiente
if max(isnan(delta)) ~= 1
 load([data_dir 'mvalold.mat']);
 temp = jacob(mvalold,theta2,args)';
 format long e
 df = 2*temp*IV*miolo*gmmresid;
end
disp('Gradient                Theta                  ');
disp('-----------------------------------------------');
disp([df theta2])
