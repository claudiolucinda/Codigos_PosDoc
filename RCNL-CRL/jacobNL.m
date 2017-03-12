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

[sij,sijg,~,~,~,numer1,denom1,numer2,denom2] = NLShareCalculation(thetaNL,Delta,Data);

% Derivative of share w.r.t. theta2

f1 = zeros(size(cdid,1),K);
%if optimalinstrument==0;

% Computing the derivative of share w.r.t. sigma

mu=mufunc(thetaNL,Data);
% part3=(sij * sij');
% part2= ((Sigseg./(1-Sigseg)).*sijg) * sij';
% part1=diag((mean(sij,2)./(1-sigsub)));
% temp_oo=dummyvar(cdid);
% 
% select2=temp_oo * temp_oo';
tam=0;
for i=1:max(Data.cdid)
    tam=tam+(size(Data.cdid(Data.cdid==i,:),1)).^2;
end

derShareDelt=spalloc(size(cdid,1),size(cdid,1),tam);

for i=1:max(Data.cdid)
    % derivative of share w.r.t. DELTA
    part3           = (sij(Data.cdid==i,:)*sij(Data.cdid==i,:)');
    part2           = ((Sigseg./(1-Sigseg)).*sijg(Data.cdid==i,:)) * sij(Data.cdid==i,:)';
    part1           = diag((sum(sij(Data.cdid==i,:),2)./(1-Sigseg)));
    temp_ooo        = sparse(dummyvar(nestid(Data.cdid==i,:)));
    Big_selector    = temp_ooo*temp_ooo';
    derShareDelt(Data.cdid==i,Data.cdid==i)    = (part1 - Big_selector.*part2 - part3)./Data.ns;
    
    % derivative os share w.r.t. sigma
    numer1t         = numer1(Data.cdid==i,:);
    numer2t         = numer2(Data.cdid==i,:);
    denom1t         = denom1(Data.cdid==i,:);
    denom2t         = denom2(Data.cdid==i,:);
    numer1t(isnan(numer1t) | isinf(numer1t))=0;
    numer2t(isnan(numer2t) | isinf(numer2t))=0;
    denom1t(isnan(denom1t) | isinf(denom1t))=0;
    denom2t(isnan(denom2t) | isinf(denom2t))=0;
    
    sgt             = numer2t./denom2t;
    
    delta_j         = (mu(Data.cdid==i,:)+(Delta(Data.cdid==i,:)*ones(1,Data.ns)));
    partA           = delta_j./((1-Sigseg).^2);
    partA(isnan(partA) | isinf(partA))  = 0;
    
    partB1          = partA.*sijg(Data.cdid==i,:);
    partB           = Sigseg.*(Big_selector*partB1);
    partB(isnan(partB) | isinf(partB))  = 0;
    
    
    partC1          = (1-Sigseg).*(partA.*sij(Data.cdid==i,:));
    partC           = ones(size(numer1t,1),1)*sum(partC1,1);  
    
    %partC1          =Big_selector*((partA.*exp(delta_j./(1-Sigseg))));
    partD           = log(denom1t);
    
    partE           = zeros(size(partD));
    for j=1:Data.ns
        temp=sum(unique((partD(:,j).*sgt(:,j)) .* ones(size(numer1t,1),1)),1);
        partE(:,j)=temp;
    end
    partE(isnan(partE) | isinf(partE))  = 0;
    
    derShareSigseg1 = sij(Data.cdid==i,:).*(partA-partB-partC-partD+partE);
    f1(Data.cdid==i,K) = mean(derShareSigseg1,2);
    
    % Derivative of share w.r.t. theta2 (RC Other than sigma)
    
    if isempty(dfull)==0;
        for m = 2:K-1
            l=m-1;
            xv            = (x2(Data.cdid==i,m)*ones(1,ns)).*vfull(Data.cdid==i,(l-1)*ns+1:l*ns);
            part1         = xv./(1-Sigseg);
            part2A        = part1.*sijg(Data.cdid==i,:);
            part2         = Sigseg.*(Big_selector*part2A);
            part2(isnan(part2) | isinf(part2))  = 0;
            part3A        = part1.*sij(Data.cdid==i,:);
            part3B        = ones(size(numer1t,1),1)*sum(part3A,1);
            part3         = (1-Sigseg).*part3B;
%             part3A        = sgt.*(Big_selector*part2A);
%             part3         = zeros(size(part3A));
%             for j=1:Data.ns
%                 temp      =sum(unique((part3A(:,j)) .* ones(size(numer1t,1),1)),1);
%                 part3(:,j)=(1-Sigseg).*temp;
%             end
            part3(isnan(part3) | isinf(part3))  = 0;
            temp          = sij(Data.cdid==i,:).*(part1-part2+part3);
            f1(Data.cdid==i,m)       = mean(temp,2);
            %clear xv temp sum1
            %f2(:,i)=mean((((p-cmg)./cmg)*ones(1,ns)).*((1./mu).*vfull(:,1:ns)+(1./(1-BigTheta*shares)).*(BigTheta*(shares.*(xv-sumxvexp)))),2);
        end
        l=1;
        
        xd            = (x2(Data.cdid==i,m)*ones(1,ns)).*dfull(Data.cdid==i,(l-1)*ns+1:l*ns);
        part1         = xd./(1-Sigseg);
        part2A        = part1.*sijg(Data.cdid==i,:);
        part2         = Sigseg.*(Big_selector*part2A);
        part2(isnan(part2) | isinf(part2))  = 0;
        part3A        = part1.*sij(Data.cdid==i,:);
        part3B        = ones(size(numer1t,1),1)*sum(part3A,1);
        part3         = (1-Sigseg).*part3B;
        %             part3A        = sgt.*(Big_selector*part2A);
        %             part3         = zeros(size(part3A));
        %             for j=1:Data.ns
        %                 temp      =sum(unique((part3A(:,j)) .* ones(size(numer1t,1),1)),1);
        %                 part3(:,j)=(1-Sigseg).*temp;
        %             end
        part3(isnan(part3) | isinf(part3))  = 0;
        temp          = sij(Data.cdid==i,:).*(part1-part2+part3);
        f1(Data.cdid==i,l)       = mean(temp,2);
    else
        for m = 2:K-1
            l=m-1;
            xv            = (x2(Data.cdid==i,m)*ones(1,ns)).*vfull(Data.cdid==i,(l-1)*ns+1:l*ns);
            part1         = xv./(1-Sigseg);
            part2A        = part1.*sijg(Data.cdid==i,:);
            part2         = Sigseg.*(Big_selector*part2A);
            part2(isnan(part2) | isinf(part2))  = 0;
            part3A        = part1.*sij(Data.cdid==i,:);
            part3B        = ones(size(numer1t,1),1)*sum(part3A,1);
            part3         = (1-Sigseg).*part3B;
%             part3A        = sgt.*(Big_selector*part2A);
%             part3         = zeros(size(part3A));
%             for j=1:Data.ns
%                 temp      =sum(unique((part3A(:,j)) .* ones(size(numer1t,1),1)),1);
%                 part3(:,j)=(1-Sigseg).*temp;
%             end

%             sgt           = numer2t./denom2t;
%             part3A        = sgt.*(Big_selector*part2A);
%             part3         = zeros(size(part3A));
%             for j=1:Data.ns
%                 temp      =sum(unique((part3A(:,j)) .* ones(size(numer1t,1),1)),1);
%                 part3(:,j)=(1-Sigseg).*temp;
%             end
            part3(isnan(part3) | isinf(part3))  = 0;
            temp          = sij(Data.cdid==i,:).*(part1-part2+part3);
            f1(Data.cdid==i,m)       = mean(temp,2);
        end
    end
end
    
    
%     
%     numer1          = exp((mu(Data.cdid==i,:)+(Delta(Data.cdid==i,:)*ones(1,Data.ns)))./(1-Sigseg));
%     denom1          = Big_selector*numer1;
%     numer2          =zeros(size(numer1));
%     for j=1:Data.ns
%         temp=sum(unique((denom1(:,j).^(1-Sigseg)) .* ones(size(numer1,1),1)),1);
%         numer2(:,j)=temp;
%     end
%     %numer2          = denom1.^(1-Sigseg);
%     denom2          = 1 + numer2;
%     
%     part1           = exp((mu(Data.cdid==i,:)+(Delta(Data.cdid==i,:)*ones(1,Data.ns)))./((1-Sigseg).^2));
%     part2A           = Big_selector*(part1.*numer1);
%     part2B           = part2A./denom1;
%     part2B(isnan(part2B)) = 0;
%     part2 = - log(denom1)-Sigseg*part2B;
%     part2(isnan(part2) | isinf(part2))  = 0;
%     
%     part3A = -log(denom1)+(1-Sigseg)*part2B;
%     part3B=zeros(size(part3A));
%     for j=1:ns
%         temp=sum(unique((denom1(:,j).^(1-Sigseg)) .* part3A(:,j)),1);
%         part3B(:,j)=temp;
%     end
%     part3                               = part3B ./ denom2;
%     part3(isnan(part3))                 = 0;
%     
    

%f = zeros(size(cdid,1),size(thetaNL,1));

parfor i = 1:max(Data.cdid)
    warning('off','MATLAB:nearlySingularMatrix');
	f2 = - derShareDelt(Data.cdid==i,Data.cdid==i)\f1(Data.cdid==i,:);
    f{i}=f2;
	
end
f=real(cat(1,f{:}));

%f=-derShareDelt\f1;
end