%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Symbolic Script for the Derivatives of share w.r.t 
% coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 char, 4 goods, 1 nest

clear; clc;

syms beta1 x11 x21 x31 x41 sigma

delta1=(beta1*x11);
delta2=(beta1*x21);
delta3=(beta1*x31);
delta4=(beta1*x41);

Denom1=(exp(delta1/(1-sigma))+exp(delta2/(1-sigma)));
Denom2=(exp(delta3/(1-sigma))+exp(delta4/(1-sigma)));


s1g1=exp(delta1/(1-sigma))/Denom1;
s2g1=exp(delta2/(1-sigma))/Denom1;


s3g2=exp(delta3/(1-sigma))/Denom2;
s4g2=exp(delta4/(1-sigma))/Denom2;

sg1=(Denom1^(1-sigma))/(1+(Denom1^(1-sigma))+(Denom2^(1-sigma)));
sg2=(Denom2^(1-sigma))/(1+(Denom1^(1-sigma))+(Denom2^(1-sigma)));

s1=s1g1*sg1;
s2=s2g1*sg1;
s3=s3g2*sg2;
s4=s4g2*sg2;
diary('C:\Users\claudiolucinda\Dropbox\Pós Doc 2016\Derivs.txt')
display('Derivative w.r.t. sigma')
ds1_dsigma=diff(s1,'sigma');

display('Not Simplified')
pretty(ds1_dsigma);
display('Simplified')
pretty(simplify(ds1_dsigma))

display('Derivative w.r.t. beta');
ds1_dbeta1=diff(s1,'beta1');
display('Not Simplified')
pretty(ds1_dbeta1)
display('Simplified')
pretty(simplify(ds1_dbeta1))

