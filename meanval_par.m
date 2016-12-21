function f = meanval_par(theta2,args)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Esta função Calcula a utilidade média
% Adaptada e otimizada por Cláudio R. Lucinda
% FGV-EESP e EAESP
% 2008/09
% Outputs:
% f - Vetor com utilidades médias
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fnames=fieldnames(args);
%for i=1:length(fnames)
%    eval([fnames{i} '=args.' fnames{i} ';']);
%end

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

theta2w = full(sparse(theti,thetj,theta2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculando aqui o mu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu=mufunc(theta2,args);
expmu = exp(mu);
norm = 1;
avgnorm = 1;


cdindex2=zeros(size(x2,1),1);
cdindex2(cdindex,1)=1;

size_proc1=round(size(cdindex,1)/2);
size_proc1=cdindex(size_proc1);
size_proc2=size(cdindex2,1)-size_proc1;
% mvalold=mvalold;
% cdid=cdid;
% s_jt=s_jt;
% ns=ns;
mvalold=mvalold;
spmd
    warning off
    codist = codistributor1d(1, [size_proc1 size_proc2]);
        expmu=codistributed(expmu,codist);
        mvalold=codistributed(mvalold,codist);
        cdindex2=codistributed(cdindex2,codist);
        cdid=codistributed(cdid,codist);
        s_jt=codistributed(s_jt,codist);

i = 0;
while norm > tol && i<2500

%while norm > tol*10^(flag*floor(i/50)) & avgnorm > 1e-3*tol*10^(flag*floor(i/50))
   mvalold2=mvalold*ones(1,ns); 
   eg = expmu.*mvalold2;
   temp = cumsum(eg); 
   sum1 = temp(cdindex2==1,:);
   sum1(2:size(sum1,1),:) = sum1(2:end,:)-sum1(1:end-1,:);
   denom1 = 1./(1+sum1);
   denom = denom1(cdid,:);
   pt_1 = eg.*denom;   
   pt_2=mean(pt_1,2);
   mval = mvalold.*s_jt./pt_2; 
   t = abs(mval-mvalold);
   norm = max(t);
   %avgnorm = mean(t);
   mvalold = mval;
   i = i + 1;
   %disp(norm)
   
end
warning on
  
end
mval=gather(mval);
mvalold=gather(mvalold2);
ii=gather(i);

disp(['Nº de Iterações - 1ª metade dos mercados:  ' num2str(ii{1})])
disp(['Nº de Iterações - 2ª metade dos mercados:  ' num2str(ii{2})])

if flag == 1 && max(isnan(mval)) < 1;
   mvalold = mval;
   oldt2 = theta2;
   save([data_dir 'mvalold.mat'],'mvalold','oldt2')
end   
f = log(mval);
