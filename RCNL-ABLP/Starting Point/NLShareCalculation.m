function [sij sijg sj sjg s0 numer1 denom1 numer2 denom2]=NLShareCalculation(thetaNL,delta,Data)

% NLShareCalculation - Nested Logit share calculation
% Syntax:   [sij sijg sj sjg s0 numer1 denom1 numer2 denom2]=NLShareCalculation(thetaNL,delta,Data)
%
% Inputs:
%    input1 - Non linear parameters
%    input2 - Data
%
% Outputs:
%    sij    = individual market shares
%    sijg   = individual conditional share of choosing good j from group g
%    sj     = market shares
%    sjg    = conditional share of choosing good j from group g 
%    s0     = outside good
%    numer1 = numerator1
%    denom1 = denominator1
%    numer2 = numerator 2
%    denom2 = denominator 2
%
% Subfunctions: none
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

Sigseg          = thetaNL(1);
RandCoeff       = thetaNL(2);

% Market share calculation
mudel           = (bsxfun(@plus,delta,xvu.*RandCoeff)) ./(1-Sigseg);
numer1          = exp(mudel);

% use multidimensional arrays
numer1          = reshape(numer1,[prods,1,nmkt,Nnodes]);
numer1          = bsxfun(@times,numer1,multidumseg);
denom1          = sum(numer1,1);

numer2          = denom1.^(1-Sigseg);
denom2          = 1 + sum(numer2,2);                         % sum across segments within a market (second dimension)


sijg              = bsxfun(@rdivide,numer1,denom1);
sijg(isnan(sijg)) = 0;

sig               = bsxfun(@rdivide,numer2,denom2);
sig(isnan(sig))   = 0;

sij               = bsxfun(@times,sijg,sig);
sij(isnan(sij))   = 0;

sj                = sum((bsxfun(@times,sij,qweightrprods)),4);
sj(isnan(sj))     = 0;

sg                = sum((bsxfun(@times,sig(:,:,:,1:Nnodes),qweightrprods)),4);             % weighted sum across draws (4th dimension)
sg(isnan(sg))     = 0;

sjg               = bsxfun(@rdivide,sj,sg);
sjg(isnan(sjg))   = 0;

s0                = 1 ./ denom2;
s0                = sum((bsxfun(@times,s0(:,:,:,1:Nnodes),qweightrprods)),4);

end