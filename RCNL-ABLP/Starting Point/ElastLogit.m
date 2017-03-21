function [ElaLogit Diversion]=ElastLogit(alpha,theta,delta,characteristic,Data)

% ElastLogit - Calculate price elasticities for Logit models
% Syntax:  [ElaLogit Diversion]=ElastLogit(alpha,theta,delta,characteristic,Data)
%
% Inputs:
%    input1 - alpha = price coefficient
%    input2 - theta = random coefficient
%    input1 - delta = mean value
%    input2 - Data
%
% Outputs:
%    prods x prods elasticity matrix
%
% Subfunctions: LogitShareCalculation
%

% Author: Laura Grigolon and Frank Verboven
% August 2012;

% Unpack
prods           = Data.prods;
nmkt            = Data.nmkt;
qweightrprods   = Data.qweightrprods;

% Market Shares
[sj, sij]       = LogitShareCalculation(theta,delta,Data);

% Jacobian
% remember that sij is prod x 1 x mkt x nodes

% derivative of shares wrt to delta
part1           = bsxfun(@times,sij,eye(prods));          % diagonal in multiple dimensions
sijtransp       = permute(sij,[2 1 3 4]);
part2           = multiprod(sij,sijtransp);
derShareDeltij  = part1 - part2;
derShareDelta   = sum((bsxfun(@times,derShareDeltij,qweightrprods)),4);
  
PriceOverShare  = multiprod(1/sj,permute(characteristic,[2 1 3]));

ElaLogit      = bsxfun(@times,(derShareDelta*alpha),PriceOverShare);

OwnDeriv        = permute((reshape(diagnd(derShareDelta),[prods,1,nmkt])),[2 1 3]);
Diversion       = - bsxfun(@times,(derShareDelta),(1/OwnDeriv));

end

