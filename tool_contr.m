function [f,img,qtiles] = tool_contr(theta2,args)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Esta função é para avaliar o desempenho - graficamente, I hope
% do "Contraction Mapping" em um determinado valor
% Cláudio R. Lucinda
% FEA/RP-USP
% 2012
% Outputs:
% f - Vetor com utilidades médias
% img - figura com a evolução da tolerância
% qtiles - quantis da distribuição de utilidades médias
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


load([data_dir 'mvalold.mat'])

tol=1e-13;
flag=1;

% if max(abs(theta2-oldt2)) < 0.01;
% % 	tol = 1e-9;
%  	flag = 0;
% else
% %   	tol = 1e-6;
%  	flag = 1;
% end

theta2w = full(sparse(theti,thetj,theta2));
%theta2w=theta2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculando aqui o mu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu=mufunc(theta2,args);
expmu = exp(mu);
norm = 1;
avgnorm = 1;


normt=[];
i = 0;
while norm > tol && i<2500

%while norm > tol*10^(flag*floor(i/50)) && avgnorm > 1e-3*tol*10^(flag*floor(i/50))
   eg = expmu.*(mvalold*ones(1,ns));
   temp = cumsum(eg); 
   sum1 = temp(cdindex,:);
   sum1(2:size(sum1,1),:) = diff(sum1);
   denom1 = 1./(1+sum1);
   denom = denom1(cdid,:);
   pt_1 = eg.*denom;   
   pt_2=mean(pt_1,2);
   mval = mvalold.*s_jt./pt_2; 
   t = abs(mval-mvalold);
   norm = max(t);
   avgnorm = mean(t);
   mvalold = mval;
   normt=[normt;norm];
   i = i + 1;
end
disp(['Nº de Iterações:  ' num2str(i)])

img=plot(normt);



if flag == 1 && max(isnan(mval)) < 1;
   mvalold = mval;
   oldt2 = theta2;
   save([data_dir 'mvalold.mat'],'mvalold','oldt2')
end

qtiles=quantile(log(mval),[.10 .25 .50 .75 .90]);

f = log(mval);
