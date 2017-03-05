function f = jacobNL(Delta,thetaNL,Data)

% jacobNL - This function is calculate the derivative of the mean value
% delta wrt the parameter using the implicit function theorem (see Nevo
% http://faculty.wcas.northwestern.edu/~ane686/supplements/Ras_guide_appendix.pdf
% for the Nested Logit
% 
% jacobNL(thetaNL,delta,Data)
%
% Inputs:
%    input1 - Non linear parameters nested logit
%    input2 - Data
%
% Outputs:
%    Derivative of Delta wrt thetaNL
%
% Subfunctions: NLShareCalculation; multinv; multiprod
%

% Author: Laura Grigolon and Frank Verboven
% August 2012;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fnames=fieldnames(Data);
for i=1:length(fnames)
    eval([fnames{i} '=Data.' fnames{i} ';']);
end

Sigseg = thetaNL(end);

K=size(thetaNL,1);

[sij,sijg,~,~,~] = NLShareCalculation(thetaNL,Delta,Data);

% Derivative of share w.r.t. theta2

f1 = zeros(size(cdid,1),K);
%if optimalinstrument==0;
    
if isempty(dfull)==0;
    for i = 2:K-1
        l=i-1;
        xv = (x2(:,i)*ones(1,ns)).*vfull(:,(l-1)*ns+1:l*ns);
        temp=xv.*((sij./(1-Sigseg))-(Sigseg/(1-Sigseg)).*sijg-sijg.^2);
        
        f1(:,i)=mean(temp,2);
        %clear xv temp sum1
        %f2(:,i)=mean((((p-cmg)./cmg)*ones(1,ns)).*((1./mu).*vfull(:,1:ns)+(1./(1-BigTheta*shares)).*(BigTheta*(shares.*(xv-sumxvexp)))),2);
    end
    clear xv temp sum1

    j=1;
    
    xd=(x2(:,j)*ones(1,ns)).*dfull(:,(j-1)*ns+1:j*ns);
    temp=xd.*((sij./(1-Sigseg))-(Sigseg/(1-Sigseg)).*sijg-sijg.^2);
        
    
    f1(:,j)=mean(temp,2)';
else
    for i = 1:K-1
        xv = (x2(:,i)*ones(1,ns)).*vfull(:,(l-1)*ns+1:l*ns);
        temp=xv.*((sij./(1-Sigseg))-(Sigseg/(1-Sigseg)).*sijg-sijg.^2);
        
        f1(:,i)=mean(temp,2);
    end
    clear xv temp sum1
end    

% Computing the derivative of share w.r.t. sigma

mu=mufunc(thetaNL,Data);
% part3=(sij * sij');
% part2= ((Sigseg./(1-Sigseg)).*sijg) * sij';
% part1=diag((mean(sij,2)./(1-sigsub)));
% temp_oo=dummyvar(cdid);
% 
% select2=temp_oo * temp_oo';

derShareDelt=sparse(size(cdid,1),size(cdid,1));

for i=1:max(Data.cdid)
    % derivative of share w.r.t. delta
    part3           = (sij(Data.cdid==i,:)*sij(Data.cdid==i,:)');
    part2           = ((Sigseg./(1-Sigseg)).*sijg(Data.cdid==i,:)) * sij(Data.cdid==i,:)';
    part1           = diag((sum(sij(Data.cdid==i,:),2)./(1-Sigseg)));
    temp_ooo        = sparse(dummyvar(nestid(Data.cdid==i,:)));
    Big_selector    = temp_ooo*temp_ooo';
    derShareDelt(Data.cdid==i,Data.cdid==i)    = (part1 - Big_selector.*part2 - part3)./ns;

    
    numer1          = exp((mu(Data.cdid==i,:)+(Delta(Data.cdid==i,:)*ones(1,Data.ns)))./(1-Sigseg));
    denom1          = Big_selector*numer1;
    numer2          =zeros(size(numer1));
    for j=1:Data.ns
        temp=sum(unique((denom1(:,j).^(1-Sigseg)) .* ones(size(numer1,1),1)),1);
        numer2(:,j)=temp;
    end
    %numer2          = denom1.^(1-Sigseg);
    denom2          = 1 + numer2;
    
    part1           = exp((mu(Data.cdid==i,:)+(Delta(Data.cdid==i,:)*ones(1,Data.ns)))./((1-Sigseg).^2));
    part2A           = Big_selector*part1;
    %part2A          = sum((part1 .* numer1),1);
    part2B          = part2A./denom1;
    part2B(isnan(part2B)) = 0;
    part2 = - log(denom1)-Sigseg*part2B;
    part2(isnan(part2) | isinf(part2))  = 0;
    
    part3A = -log(denom1)+(1-Sigseg)*part2B;
    part3B=zeros(size(part3A));
    for j=1:ns
        temp=sum(unique((denom1(:,j).^(1-Sigseg)) .* part3A(:,j)),1);
        part3B(:,j)=temp;
    end
    part3                               = part3B ./ denom2;
    part3(isnan(part3))                 = 0;
    derShareSigseg1                     = sij(Data.cdid==i,:).*(part1+part2-part3);
    f1(Data.cdid==i,K) = mean(derShareSigseg1,2);
    
end

f = zeros(size(cdid,1),size(thetaNL,1));

for i = 1:size(Data.cdid,1)
	f(Data.cdid==i,:) = - derShareDelt(Data.cdid==i,Data.cdid==i)\f1(Data.cdid==i,:);
	
end


%f=-derShareDelt\f1;
end