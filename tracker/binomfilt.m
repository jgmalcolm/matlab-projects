function [h] = binomfilt(order1,order2)

% BINOMFILT Binomial Kernel
%   BIMONFILT(order1) constructs binomial kernel of
%   order n = order1. 
%   BIMONFILT(order1,order2) constructs a separable
%   2-D binomial kernel as a tensor (outer) product
%   of two 1-D kernels having the orders n = order1
%   and m = order2, respectively.
%
% written by Oleg Michailovich, July 2004

h=zeros(order1+1,1);
for k=0:order1
    h(k+1)=nchoosek(order1,k);
end
h=h/(2^order1);

if (nargin>1)
    g=zeros(order2+1,1);
    for k=0:order2
        g(k+1)=nchoosek(order2,k);
    end
    g=g/(2^order2);
    h=h*g';
end