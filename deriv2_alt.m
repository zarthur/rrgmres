function [A,b,x] = deriv2(n,example)
%DERIV2 Test problem: computation of the second derivative.
%
% [A,b,x] = deriv2_alt(n,example)
%
% This is a mildly ill-posed problem.  It is a discretization of a
% first kind Fredholm integral equation whose kernel K is the
% Green's function for the second derivative:
%    K(s,t) = | s(t-1)  ,  s <  t .
%             | t(s-1)  ,  s >= t
% Both integration intervals are [0,1], and as right-hand side g
% and correspond solution f one can choose between the following:
%    example = 1 : g(s) = (s^3 - s)/6          ,  f(t) = t
%    example = 2 : g(s) = exp(s) + (1-e)s - 1  ,  f(t) = exp(t)
%    example = 3 : g(s) = | (4s^3 - 3s)/24               ,  s <  0.5
%                         | (-4s^3 + 12s^2 - 9s + 1)/24  ,  s >= 0.5
%                  f(t) = | t    ,  t <  0.5
%                         | 1-t  ,  t >= 0.5

% References.  The first two examples are from L. M. Delves & J. L.
% Mohamed, "Computational Methods for Integral Equations", Cambridge
% University Press, 1985; p. 310.  The third example is from A. K.
% Louis & P. Maass, "A mollifier method for linear operator equations
% of the first kind", Inverse Problems 6 (1990), 427-440.

% Discretized by the trapezoid rule.

% Initialization.
if (nargin==1), example = 1; end
h = 1/(n-1); A = zeros(n,n); h2 = h^2;

% Compute the matrix A.
for i = 1:n
  A(i,i) = h2*(i^2 - (n+1)*i + n);
  for j = (i+1):n
    A(i,j) = h2*(i-1)*(j-n);
    A(j,i) = A(i,j);
  end
end
for i = 1:n
    A(i,1) = A(i,1)/2;
    A(i,n) = A(i,n)/2;
end

% Compute the right-hand side vector b.
if (nargout>1)
    b = zeros(n,1);
    if (example==1)
        for i=1:n
            s_i = (i-1)*h;
            b(i) = (s_i^3 - s_i)/6;
        end
    elseif (example==2)
        for i=1:n
            s_i = (i-1)*h;
            b(i) = exp(s_i) + (1-exp(1))*s_i - 1;
        end
    else
        error('Illegal value of example')
    end
end

% Compute the solution vector x.
if (nargout==3)
    x = zeros(n,1);
    if (example==1)
        for i=1:n
            s_i = (i-1)*h;
            x(i) = s_i;
        end
    elseif (example==2)
        for i=1:n
            s_i = (i-1)*h;
            x(i) = exp(s_i);
        end
    else
        error('Illegal value of example')
    end
end

