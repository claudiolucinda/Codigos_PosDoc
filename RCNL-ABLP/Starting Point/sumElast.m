function [MeanOwnElast OneMeanOwnElast MeanCrossESeg OneMeanCrossESeg MeanCrossEDiff OneMeanCrossEDiff MeanDiversion MeanDiversionDiff OneMeanDiversion OneMeanDiversionDiff] = sumElast(Elasticity,Diversion,Data)

% sumElast - Summarize Elasticities by segment
% Syntax:  [MeanOwnElast OneMeanOwnElast MeanCrossESeg OneMeanCrossESeg MeanCrossEDiff OneMeanCrossEDiff MeanDiversion MeanDiversionDiff OneMeanDiversion OneMeanDiversionDiff] = sumElast(Elasticity,Diversion,Data)
%
% Inputs:
%    Elasticity
%    input2 - Data
%
% Outputs:
%    Mean Value Own Elasticity, Mean Value Cross Elasticity wrt same
%    segment; Mean Value Cross Elasticity wrt different segment
%
% Subfunctions: diagnd
%

% Author: Laura Grigolon and Frank Verboven
% August 2012;

prods           = Data.prods;
nmkt            = Data.nmkt;
multidumseg     = Data.multidumseg;

multidumElast     = multidumseg(:,:,:,1);    % retain only 3 dimensions
OwnElast          = bsxfun(@times,reshape(diagnd(Elasticity),[prods,1,nmkt]),multidumElast);
MeanOwnElast      = sum(OwnElast)./sum(OwnElast~=0);

OneMeanOwnElast   = mean(sum(OwnElast,2));    

multidumsegtrans  = permute(multidumElast,[2 1 3]);
SelEla            = multiprod(multidumElast,multidumsegtrans);
SelEla1           = bsxfun(@times,multidumElast(:,1,:),SelEla);
SelEla2           = bsxfun(@times,multidumElast(:,2,:),SelEla);
CrossESameSeg1    = bsxfun(@times,(bsxfun(@times,SelEla1,Elasticity)),1-eye(prods)); % 1-eye is just to get rid of the diagonal
CrossESameSeg2    = bsxfun(@times,(bsxfun(@times,SelEla2,Elasticity)),1-eye(prods));
MeanCrossESeg     = [sum(sum(CrossESameSeg1))./sum(sum(CrossESameSeg1~=0)),...
                     sum(sum(CrossESameSeg2))./sum(sum(CrossESameSeg2~=0))];

                 
OneMeanCrossESeg  = (sum(sum(CrossESameSeg1)) + sum(sum(CrossESameSeg2)))./...
                    (sum(sum(CrossESameSeg1~=0)) + sum(sum(CrossESameSeg2~=0)));

DiversionSeg1     = bsxfun(@times,(bsxfun(@times,SelEla1,Diversion)),1-eye(prods)); % 1-eye is just to get rid of the diagonal
DiversionSeg2     = bsxfun(@times,(bsxfun(@times,SelEla2,Diversion)),1-eye(prods)); % 1-eye is just to get rid of the diagonal

MeanDiversion     = [sum(sum(DiversionSeg1))./sum(sum(DiversionSeg1~=0)),...
                     sum(sum(DiversionSeg2))./sum(sum(DiversionSeg2~=0))];
                 
OneMeanDiversion  = (sum(sum(DiversionSeg1)) + sum(sum(DiversionSeg2)))./...
                     (sum(sum(DiversionSeg1~=0)) + sum(sum(DiversionSeg2~=0)));

MirrorSelEla      = 1-SelEla;
SelEla1           = bsxfun(@times,multidumElast(:,1,:),MirrorSelEla);
SelEla2           = bsxfun(@times,multidumElast(:,2,:),MirrorSelEla);
CrossEDiffSeg1    = bsxfun(@times,(bsxfun(@times,SelEla1,Elasticity)),1-eye(prods)); % 1-eye is just to get rid of the diagonal
CrossEDiffSeg2    = bsxfun(@times,(bsxfun(@times,SelEla2,Elasticity)),1-eye(prods));
MeanCrossEDiff    = [sum(sum(CrossEDiffSeg1))./sum(sum(CrossEDiffSeg1~=0)),...
                     sum(sum(CrossEDiffSeg2))./sum(sum(CrossEDiffSeg2~=0))];
                 
OneMeanCrossEDiff  = (sum(sum(CrossEDiffSeg1)) + sum(sum(CrossEDiffSeg2)))./...
                     (sum(sum(CrossEDiffSeg1~=0)) + sum(sum(CrossEDiffSeg2~=0)));

DiversionDiffSeg1  = bsxfun(@times,SelEla1,Diversion); % 1-eye is just to get rid of the diagonal
DiversionDiffSeg2  = bsxfun(@times,SelEla2,Diversion); % 1-eye is just to get rid of the diagonal

MeanDiversionDiff  =  [sum(sum(DiversionDiffSeg1))./sum(sum(DiversionDiffSeg1~=0)),...
                      sum(sum(DiversionDiffSeg2))./sum(sum(DiversionDiffSeg2~=0))];
                  
OneMeanDiversionDiff  = (sum(sum(DiversionDiffSeg1)) + sum(sum(DiversionDiffSeg2)))./...
                        (sum(sum(DiversionSeg1~=0)) + sum(sum(DiversionSeg2~=0)));

end