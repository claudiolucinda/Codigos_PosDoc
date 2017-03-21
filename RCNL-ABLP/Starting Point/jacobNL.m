function f = jacobNL(thetaNL,delta,Data)

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

% Unpack
prods           = Data.prods;
nmkt            = Data.nmkt;
Nnodes          = Data.Nnodes;
qweightrprods   = Data.qweightrprods;
multidumseg     = Data.multidumseg;
xvu             = Data.xvu;
multixvu        = Data.multixvu;

% Market Shares
[sij,sijg,~,~,~, numer1,denom1,~,denom2] = NLShareCalculation(thetaNL,delta,Data);
Sigseg          = thetaNL(1);
RandCoeff       = thetaNL(2);
% Derivative wrt delta

part1          = sum((1./(1-Sigseg) * sij),2);
part1          = bsxfun(@times,part1,eye(prods));          % diagonal in multiple dimensions

sijtransp      = permute(sij,[2 1 3 4]);
part2          = (Sigseg/(1-Sigseg))*multiprod(sijg,sijtransp);

part3          = multiprod(sum(sij,2),sum(sijtransp,1));

derShareDeltij = (part1 - part2 - part3);
derShareDelta  = sum((bsxfun(@times,derShareDeltij,qweightrprods)),4);

% Derivative wrt the random coefficients

part1               = multixvu ./(1-Sigseg);
part2A              = sum((bsxfun(@times,numer1,part1)),1);
part2               = Sigseg * (part2A ./denom1);
part2(isnan(part2)) = 0;

part3A              = (1-Sigseg) * ((denom1 .^ (-Sigseg)) .* (part2A));
part3B              = sum(part3A,2);
part3               = part3B ./ denom2;
part3(isnan(part3)) = 0;

derShareRC1    = bsxfun(@minus,part1,part2);
derShareRC2    = sij .* bsxfun(@minus,derShareRC1,part3);
derShareRC     = sum((bsxfun(@times,derShareRC2,qweightrprods)),4);

% Derivative wrt the sigma segment

part1                               = (bsxfun(@plus,delta,(xvu.*RandCoeff))) ./ ((1-Sigseg)^2);
part1                               = reshape(part1,[prods,1,nmkt,Nnodes]);
part1                               = bsxfun(@times,part1,multidumseg);

part2A                              = sum((part1 .* numer1),1);
part2B                              = bsxfun(@rdivide,part2A,denom1);
part2B(isnan(part2B))               = 0;
part2                               = bsxfun(@minus,(- log(denom1)),(Sigseg *part2B));
part2(isnan(part2) | isinf(part2))  = 0;

part3A                              = bsxfun(@plus,(- log(denom1)),((1-Sigseg) *part2B));
part3B                              = sum(((denom1.^(1-Sigseg)) .* part3A),2);
part3                               = part3B ./ denom2;
part3(isnan(part3))                 = 0;

derShareSigseg1                     = bsxfun(@minus,part2,part3);
derShareSigseg2                     = sij .* bsxfun(@plus,part1,derShareSigseg1);
derShareSigseg                      = sum((bsxfun(@times,derShareSigseg2,qweightrprods)),4);


derShareTheta                       = [sum(derShareSigseg,2) , sum(derShareRC,2)];
derDeltaTheta                       = multiprod((-multinv(derShareDelta)),derShareTheta);
derDeltaTheta                       = [reshape(derDeltaTheta(:,1,:,:),[Data.nobs 1]) , reshape(derDeltaTheta(:,2,:,:),[Data.nobs 1])]; % 2D market shares

f = derDeltaTheta;