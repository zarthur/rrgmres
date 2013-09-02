function [A,b,x] = baart_alt(n) 
%BAART Test problem: Fredholm integral equation of the first kind. 
% 
% [A,b,x] = baart(n) 
% 
% Discretization of a first-kind Fredholm integral equation with 
% kernel K and right-hand side g given by 
%    K(s,t) = exp(s*cos(t)) ,  g(s) = 2*sinh(s)/s , 
% and with integration intervals  s in [0,pi/2] ,  t in [0,pi] . 
% The solution is given by 
%    f(t) = sin(t) . 
% 
% The order n must be even. 
  
% Reference: M. L. Baart, "The use of auto-correlation for pseudo- 
% rank determination in noisy ill-conditioned linear least-squares 
% problems", IMA J. Numer. Anal. 2 (1982), 241-247. 
 
% Discretized by trapezoid rule;
 
 
% Check input. 
if (n < 1), error('The order n must be positive'), end 
 
% Generate the matrix. 
hs = pi/(2*(n-1)); ht = pi/(n-1); A = zeros(n,n);

for i = 1:n
    hsi = (i-1)*hs;
    for j = 1:n
        htj = (j-1)*ht;
        A(i,j) = ht*exp(hsi*cos(htj));
    end
end 
for i = 1:n
    A(i,1)=A(i,1)/2;
    A(i,n)=A(i,n)/2;
end
% Generate the right-hand side. 
if (nargout>1)
    b = zeros(n,1);
    b(1)=2;
    for i = 2:n
        hsi = (i-1)*hs;
        b(i) = 2*sinh(hsi)/hsi;
    end
end
 
% Generate the solution. 
if (nargout==3)
    x = zeros(n,1);
    for j = 1:n
        htj = (j-1)*ht;
        x(j) = sin(htj);
    end
end 