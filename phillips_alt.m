function [A,b,xx] = phillips_alt(n);
%PHILLIPS "famous" test problem.
%
% [A,b,x] = phillips_alt(n)
%
% Discretization of the `famous' first-kind Fredholm integral
% equation deviced by D. L. Phillips.  Define the function
%    phi(x) = | 1 + cos(x*pi/3) ,  |x| <  3 .
%             | 0               ,  |x| >= 3
% Then the kernel K, the solution f, and the right-hand side
% g are given by:
%    K(s,t) = phi(s-t) ,
%    f(t)   = phi(t) ,
%    g(s)   = (6-|s|)*(1+.5*cos(s*pi/3)) + 9/(2*pi)*sin(|s|*pi/3) .
% Both integration intervals are [-6,6].
%
% Reference: D. L. Phillips, "A technique for the numerical solution
% of certain integral equations of the first kind", J. ACM 9
% (1962), 84-97.
%
% Discretized by trapezoid rule.
%
% Check input.
if (n < 2), error('The order n must be at least 2'), end
%
% Compute the matrix A.
h = 12/(n-1); r1 = zeros(1,n);
if (rem((n-1),4) == 0)
    for i = 1:((n-1)/4)
        r1(i) = h*(1 + cos((i-1)*h*pi/3));
    end
else
    for i = 1:(fix((n-1)/4)+1)
        r1(i) = h*(1 + cos((i-1)*h*pi/3));
    end
end
A = toeplitz(r1);
for i = 1:(fix((n-1)/4)+1)
    A(i,1) = A(i,1)/2;
    A(n-i+1,n) = A(n-i+1,n)/2;
end

% Compute the right-hand side b.
if (nargout>1),
  b = zeros(n,1);
  for i = 1:n
      s = -6 + (i - 1)*h;
      b(i) = (6 - abs(s))*(1+0.5*cos(s*pi/3)) + 9/(2*pi)*sin(abs(s)*pi/3);
  end
end

% Compute the solution x.
if (nargout==3),
  xx = zeros(n,1);
  i1 = fix((n-1)/4)+2;
  i2 = fix(3*(n-1)/4)+1;
  for i = i1:i2
      s = -6 + (i - 1)*h;
      xx(i)= 1 + cos(s*pi/3);
  end
end
%cond(A)