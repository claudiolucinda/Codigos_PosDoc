function [ElaNLogit Diversion]=ElastNestedLogit(alpha,thetaNL,delta,characteristic,Data)

% ElastNestedLogit - Calculate price elasticities for Nested Logit models
% Syntax:  f=ElastLogit(alpha,theta,delta,Data)
%
% Inputs:
%    input1 - alpha     = price coefficient
%    input2 - thetaNL   = random coefficient
%    input1 - delta     = mean value
%    input2 - Data
%
% Outputs:
%    prods x prods elasticity matrix
%
% Subfunctions: NLShareCalculation 
%
% Author: Laura Grigolon and Frank Verboven
% August 2012;

% Unpack
prods             = Data.prods;
nmkt             = Data.nmkt;
qweightrprods     = Data.qweightrprods;

% Market Shares
[sij,sijg,sj,~,~, ~,~,~,~] = NLShareCalculation(thetaNL,delta,Data);
Sigseg          = thetaNL(1);

% Derivative wrt delta
part1            = sum((1./(1-Sigseg) * sij),2);
part1            = bsxfun(@times,part1,eye(prods));          % diagonal in multiple dimensions

sijtransp        = permute(sij,[2 1 3 4]);
part2            = (Sigseg/(1-Sigseg))*multiprod(sijg,sijtransp);

part3            = multiprod(sum(sij,2),sum(sijtransp,1));

derShareDeltij   = (part1 - part2 - part3);
derShareDelta    = sum((bsxfun(@times,derShareDeltij,qweightrprods)),4);

PriceOverShare   = multiprod((1/(sum(sj,2))),permute(characteristic,[2 1 3]));

ElaNLogit = bsxfun(@times,(derShareDelta*alpha),PriceOverShare);

OwnDeriv        = permute((reshape(diagnd(derShareDelta),[prods,1,nmkt])),[2 1 3]);
Diversion       = bsxfun(@times,(derShareDelta),(1/OwnDeriv));

end
