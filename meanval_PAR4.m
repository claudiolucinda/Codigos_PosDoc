function f = meanval_PAR4(theta2,args)
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

%global args

fnames=fieldnames(args);
for i=1:length(fnames)
    eval([fnames{i} '=args.' fnames{i} ';']);
end

load([data_dir 'mvalold.mat'])

tol=1e-13;
flag=1;

theta2w = full(sparse(theti,thetj,theta2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculando aqui o mu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu=mufunc(theta2,args);
expmu = exp(mu);

matlabpool size;
nworkers=ans;
temp3=1;
seletor=1:size(s_jt,1);
for t=1:nworkers-1
    temp1=floor(t*size(cdindex,1)/nworkers);
    temp2=cdindex(temp1);
	mktsubs{t} =  seletor>=temp3 & seletor<=temp2;
 	mkts_jt{t} = s_jt(mktsubs{t},:);
	mktmvalold{t} = mvalold(mktsubs{t},:);
	mktexpmu{t} = expmu(mktsubs{t},:);
    mktcdid{t} = cdid(mktsubs{t},:);
    temp3=cdindex(temp1)+1;
    nsM{t} = ns;
    
end
t=nworkers;
temp2=cdindex(end);
mktsubs{t} =  seletor>=temp3 & seletor<=temp2;
mkts_jt{t} = s_jt(mktsubs{t},:);
mktmvalold{t} = mvalold(mktsubs{t},:);
mktexpmu{t} = expmu(mktsubs{t},:);
mktcdid{t} = cdid(mktsubs{t},:);
nsM{t} = ns;

%eval(['amost' num2str(nworkers) '=zeros(size(s_jt));']);
%eval(['amost' num2str(nworkers) '(temp3:temp2,:)=1;']);
%mktsubs{nworkers} =  find(cdid>temp3 & cdid<temp2);

warning off
parfor t=1:nworkers,
    %address=eval(['find(amost' num2str(t) '==1);']);
    share=mkts_jt{t};
    mvalmval=mktmvalold{t};
    expmu_m=mktexpmu{t};
    cdid_m=mktcdid{t};
    cdid_m=(cdid_m-min(cdid_m))+1;
    temp_oo=dummyvar(cdid_m);
    ns_m=nsM{t};
    
        
%     share=mkts_jt{t};
%     mvalmval=mktmvalold{t};
%     expmu_m=mktexpmu{t};
%     cdid_m=mktcdid{t};
%     cdid_m=(cdid_m-min(cdid_m))+1;
 
    % compute exp(mu) once
    %emu=exp(mkt.x*(draws.nu.*repmat(theta2',[ns 1]))');
    % iterate contraction mapping to convergence
    i2=0;
    norm=1;
    while norm > tol && i2<2500
        i2 = i2+1; 
        eg = expmu_m.*(mvalmval*ones(1,ns_m));
        sumexp=temp_oo'*eg;
        denom = 1./(1+sumexp);
        sum1=temp_oo*denom;
        
        pt_1 = eg.*sum1;   
        pt_2 = mean(pt_1,2);
        mval = mvalmval.*share./pt_2;
        difdif = abs(mval-mvalmval);
        norm = max(difdif);
        mvalmval = mval;
%        disp(['norma - parte' num2str(t) ' iter ' num2str(i2) '=' num2str(norm)]);
        
    end
    mval_par{t}=num2cell(mval);
    disp(['N�mero de Itera��es na parte ' num2str(t) ' dos mercados: ' num2str(i2)]);
end
warning on
    
f=zeros(size(s_jt));
for t=1:nworkers
   f(mktsubs{t},:)=cell2mat(mval_par{t});
end
f=log(f);

if flag == 1 && max(isnan(f)) < 1;
   mvalold = f;
   %oldt2 = theta2;
   save([data_dir 'mvalold.mat'],'mvalold')
   %save([data_dir 'mvalold.mat'],'mvalold','oldt2')
end   




% 
% spmd
%     codist = codistributor1d(1, [size_proc1 size_proc2]);
%     codist2 = codistributor1d(1,[1 1]);
%     i1=codistributed(i0,codist2);
%     i2=getLocalPart(i1);
%     expmu=codistributed(expmu,codist);
%     expmu2=getLocalPart(expmu);
%     mvalold=codistributed(mvalold,codist);
%     mvalold2=getLocalPart(mvalold);
%     cdindex2=codistributed(cdindex2,codist);
%     cdindex22=getLocalPart(cdindex2);
%     cdid=codistributed(cdid,codist);
%     cdid2=getLocalPart(cdid);
%     s_jt=codistributed(s_jt,codist);
%     s_jt2=getLocalPart(s_jt);
%         
%        
% 
% 
% while norm > tol && i2<2500
% 
%    mvalold22=mvalold2*ones(1,ns); 
%    eg = expmu2.*mvalold22;
%    temp_oo=dummyvar(cdid2);
%    sumshares=temp_oo'*eg;
%    sum1=temp_oo*sumshares;
%    denom1 = 1./(1+sum1);
%    denom = denom1(cdid2,:);
%    pt_1 = eg.*denom;   
%    pt_2=mean(pt_1,2);
%    mval = mvalold2.*s_jt2./pt_2; 
%    t = abs(mval-mvalold2);
%    norm = max(t);
%    mvalold2 = mval;
%    i2 = i2 + 1;
%    
%    
%    
%   
% end
% %warning on
% end
% mval=gather(mval);
% mval=mval{:};
% mvalold=gather(mvalold);
% ii=gather(i2);
% 
% disp(['N� de Itera��es - 1� metade dos mercados:  ' num2str(ii{1})])
% disp(['N� de Itera��es - 2� metade dos mercados:  ' num2str(ii{2})])
% if flag == 1 && max(isnan(mval)) < 1;
%    mvalold = mval;
%    oldt2 = theta2;
%    save([data_dir 'mvalold.mat'],'mvalold','oldt2')
% end   
% f = log(mval);
