function [deltas,shares,sharecond] = testerLOOP(thetaNL,Data)

% DeltaNL - Contraction Mapping for Random Coefficient Nested Logit
% Works like Logit apart for the weighting of the difference between
% calculated and observed shares
% Syntax:  f = deltaNL(thetaNL,Data)
%
% Inputs:
%    input1 - Non linear parameters
%    input2 - Data
%
% Outputs:
%    Mean Value Delta
%
% MAT-files required: mvaloldNL
% Subfunctions: NLShareCalculation
%

% Author: Laura Grigolon and Frank Verboven
% August 2012;


fnames=fieldnames(Data);
for i=1:length(fnames)
    eval([fnames{i} '=Data.' fnames{i} ';']);
end

load([data_dir 'mvaloldNL.mat'],'delta0NL');
if max(isnan(delta0NL)|isinf(delta0NL))  == 1
    deltastart           =  zeros(size(delta0NL));
else
    deltastart           = delta0NL;
end

tol=1e-13;

%logobsshare              = log(sj);


warning off
% Fazendo Mercado Por Mercado Paralelizando
mu = mufunc(thetaNL,Data);
Sigseg = thetaNL(end);
disp('Contraction Mapping - Time:');
tic
i=1;
deltas=[];

%parfor i=1:max(Data.cdid)
    deltastart=delta0NL(Data.cdid==i,:);
    deltas=[deltas deltastart];
    share=sj(Data.cdid==i,:);
    i2=0;
    norm=1;
    shares=share;
    sharecond=[];
    while norm > tol && i2<100
        i2 = i2+1;
        
        expDelta=deltastart*ones(1,Data.ns);
        % Market share calculation
        
        numer1          = exp((mu(Data.cdid==i,:)+expDelta)./(1-Sigseg));
        
        temp_ooo=sparse(dummyvar(Data.nestid(Data.cdid==i,:)));
        matrix_sel2=sparse(temp_ooo*temp_ooo');
        denom1=matrix_sel2*(numer1);
        numer2=denom1.^(1-Sigseg);

        denom2=zeros(size(denom1));

        % Versão sem paralelização
        for j=1:Data.ns
            temp=sum(unique(numer2(:,j)),1);
            denom2(:,j)=temp*ones(size(denom2(:,j)));
        end
        
        denom2=1+denom2;

        sijg=numer1./denom1;
        sijg(isnan(sijg)) = 0;

        sig=numer2./denom2;
        sig(isnan(sig))   = 0;

        sij=sig.*sijg;
        sij(isnan(sij)) = 0;

        sh=mean(sij,2);
        sh(isnan(sh)) = 0;
        
        shares=[shares sh];
        sharecond=[sharecond mean(sijg,2)];
        
        %delta1  = deltastart.*((share./sh).^(1-Sigseg));
        delta1  = deltastart+(1-Sigseg)*(log(share)-log(sh));
        
        norm = max(abs(delta1-deltastart));
        disp('Max. Tol.');
        disp(num2str(max(abs(delta1-deltastart))))
        disp('Mean. Tol.');
        disp(num2str(mean(abs(delta1-deltastart),1)))
        deltastart  = delta1;
        deltas=[deltas deltastart];

        
    end

