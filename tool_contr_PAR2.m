function [f,iters,norms] = tool_contr_PAR2(theta2,args)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Esta função dá mais informações sobre os contraction mappings
% Adaptada e otimizada por Cláudio R. Lucinda
% Dicas importantes de um Código de Conlon (2011)
% FEA-RP/USP
% 2013
% Outputs:
% f - Vetor com utilidades médias
% iters - vetor com o número de iterações nos mercadps
% mval - Valor do abs(mval-mvalold) no final
% Inputs:
% theta2 - Vetor de Parametros (o theta2w "amassado")
% args - estrutura com os seguinte campos
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

%matlabpool size;
nworkers=max(cdid);
%temp3=1;
%seletor=1:size(s_jt,1);
for t=1:nworkers
    %temp1=floor(t*size(cdindex,1)/nworkers);
    %temp2=cdindex(temp1)-1;
    %eval(['amost' num2str(t) '=zeros(size(s_jt));']);
    %eval(['amost' num2str(t) '(temp3:temp2,:)=1;']);
	mktsubs{t} =  cdid==t;
 	mkts_jt{t} = s_jt(mktsubs{t},:);
	mktmvalold{t} = mvalold(mktsubs{t},:);
	mktexpmu{t} = expmu(mktsubs{t},:);
    mktcdid{t} = cdid(mktsubs{t},:);
%    temp3=cdindex(temp1);
end

% t=nworkers;
% temp2=cdindex(end);
% mktsubs{t} =  seletor>temp3 & seletor<temp2;
% mkts_jt{t} = s_jt(mktsubs{t},:);
% mktmvalold{t} = mvalold(mktsubs{t},:);
% mktexpmu{t} = expmu(mktsubs{t},:);
% mktcdid{t} = cdid(mktsubs{t},:);

%eval(['amost' num2str(nworkers) '=zeros(size(s_jt));']);
%eval(['amost' num2str(nworkers) '(temp3:temp2,:)=1;']);
%mktsubs{nworkers} =  find(cdid>temp3 & cdid<temp2);

warning off
for t=1:nworkers,
    %address=eval(['find(amost' num2str(t) '==1);']);
    share=mkts_jt{t};
    mvalmval=mktmvalold{t};
    expmu_m=mktexpmu{t};
    cdid_m=mktcdid{t};
    cdid_m=(cdid_m-min(cdid_m))+1;
    temp_oo=dummyvar(cdid_m);
        
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
        eg = expmu_m.*(mvalmval*ones(1,ns));
        sumexp=temp_oo'*eg;
        denom = 1./(1+sumexp);
        sum1=temp_oo*denom;
        pt_1 = eg.*sum1;   
        pt_2 = mean(pt_1,2);
        mval = mvalmval.*share./pt_2;
        difdif = abs(mval-mvalmval);
        norm = max(difdif);
        mvalmval = mval;
        i2 = i2+1; 
        disp(['norma - parte' num2str(t) ' iter ' num2str(i2) '=' num2str(norm)]);
        
    end
    it_it{t}=i2;
    norm_norm{t}=norm;
    mval_par{t}=num2cell(mval);
    disp(['Número de Iterações na parte ' num2str(t) ' dos mercados: ' num2str(i2)]);
    disp(['delta médio: ' num2str(mean(mval,1))]);
        
end
warning on

norms=zeros(nworkers,1);
iters=zeros(nworkers,1);
f=zeros(size(s_jt));
for t=1:nworkers
    norms(t,1)=norm_norm{t};
    iters(t,1)=it_it{t};
   f(mktsubs{t},:)=cell2mat(mval_par{t});
end
clear temp

if flag == 1 && max(isnan(f)) < 1;
   mvalold = f;
   oldt2 = theta2;
   save([data_dir 'mvalold.mat'],'mvalold','oldt2')
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
% disp(['Nº de Iterações - 1ª metade dos mercados:  ' num2str(ii{1})])
% disp(['Nº de Iterações - 2ª metade dos mercados:  ' num2str(ii{2})])
% if flag == 1 && max(isnan(mval)) < 1;
%    mvalold = mval;
%    oldt2 = theta2;
%    save([data_dir 'mvalold.mat'],'mvalold','oldt2')
% end   
% f = log(mval);
