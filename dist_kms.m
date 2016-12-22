%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to generate Random Draws from the empirical distribution of the 
% Mileage - São Paulo
% Importance sampling approach
% Copyright 2016 Claudio R. Lucinda
% FEA-RP/USP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [results]=dist_kms(ndraws,nlines)
%clear; clc;

%ndraws=1e6;
%nlines=2;
probs=[100.00	99.25	98.50	97.75	97.00	94.60	92.20	89.80	...
87.40	85.00	79.00	73.00	67.00	61.00	55.00	50.00	45.00	40.00	...
35.00	30.00	27.60	25.20	22.80	20.40	18.00	16.40	14.80	13.20	...
11.60	10.00	9.00	8.00	7.00	6.00	5.00	4.20	3.40	2.60	...
1.80	1.00]';

kms=@(x) (11266+779.66.*x-49.566.*(x.^2)+0.6716.*(x.^3));
cum_kms=cumsum(kms(1:1:40));

%cum_probs=100-probs;
f_probs=(probs(1:end-1)-probs(2:end))./100;
f_probs=[f_probs;0.01];
q_i=randi(40,ndraws,1);
shift=1/40;
%p_i=f_probs(q_i);

w_i=f_probs./shift;

%is_dist=cum_kms(q_i).*(p_i./shift);
%y=randsample(40,ndraws,true,w_i);
[~, K] = histc(rand(ndraws,nlines),cumsum([0;w_i(:)./sum(w_i)]));
results=cum_kms(K)';
end