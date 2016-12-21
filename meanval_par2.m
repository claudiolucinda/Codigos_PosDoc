function f = meanval_par2(theta2,args)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Esta fun��o Calcula a utilidade m�dia
% Adaptada e otimizada por Cl�udio R. Lucinda
% FGV-EESP e EAESP
% 2008/09
% Outputs:
% f - Vetor com utilidades m�dias
% df - Gradiente da Fun��o - Para otimizadores que precisam
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
%   vfull - parte n�o observ�vel j� expandida
%   dfull - parte draws demografia j� expandida
%   cdindex - vetor que quebra os mercados
%   cdid - vetor que diz a que mercado pertence cada observa��o
%   data_dir - string com o diret�rio aonde est�o os dados
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

cdindex2=zeros(size(x2,1),1);
cdindex2(cdindex,1)=1;

size_proc1=round(size(cdindex,1)/2);
size_proc1=cdindex(size_proc1);
size_proc2=size(cdindex2,1)-size_proc1;
    
mvalold=mvalold;
i0=[0;0];

spmd
    codist = codistributor1d(1, [size_proc1 size_proc2]);
    codist2 = codistributor1d(1,[1 1]);
    i1=codistributed(i0,codist2);
    i2=getLocalPart(i1);
    expmu=codistributed(expmu,codist);
    expmu2=getLocalPart(expmu);
    mvalold=codistributed(mvalold,codist);
    mvalold2=getLocalPart(mvalold);
    cdindex2=codistributed(cdindex2,codist);
    cdindex22=getLocalPart(cdindex2);
    cdid=codistributed(cdid,codist);
    cdid2=getLocalPart(cdid);
    s_jt=codistributed(s_jt,codist);
    s_jt2=getLocalPart(s_jt);
        
       


while norm > tol && i2<2500

   mvalold22=mvalold2*ones(1,ns); 
   eg = expmu2.*mvalold22;
   temp_oo=dummyvar(cdid2);
   sumshares=temp_oo'*eg;
   sum1=temp_oo*sumshares;
   denom1 = 1./(1+sum1);
   denom = denom1(cdid2,:);
   pt_1 = eg.*denom;   
   pt_2=mean(pt_1,2);
   mval = mvalold2.*s_jt2./pt_2; 
   t = abs(mval-mvalold2);
   norm = max(t);
   mvalold2 = mval;
   i2 = i2 + 1;
   
   
   
  
end
%warning on
end
mval=gather(mval);
mval=mval{:};
mvalold=gather(mvalold);
ii=gather(i2);

disp(['N� de Itera��es - 1� metade dos mercados:  ' num2str(ii{1})])
disp(['N� de Itera��es - 2� metade dos mercados:  ' num2str(ii{2})])
if flag == 1 && max(isnan(mval)) < 1;
   mvalold = mval;
   oldt2 = theta2;
   save([data_dir 'mvalold.mat'],'mvalold','oldt2')
end   
f = log(mval);
