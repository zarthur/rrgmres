function [A,b,x] = shaw_alt(n)
%SHAW Test problem: one-dimensional image restoration model.
%
% [A,b,x] = shaw(n)
%
% Discretization of a first kind Fredholm integral equation with
% [-pi/2,pi/2] as both integration intervals.  The kernel K and
% the solution f, which are given by
%    K(s,t) = (cos(s) + cos(t))*(sin(u)/u)^2
%    u = pi*(sin(s) + sin(t))
%    f(t) = a1*exp(-c1*(t - t1)^2) + a2*exp(-c2*(t - t2)^2) ,
% are discretized by simple quadrature to produce A and x.
% Then the discrete right-hand b side is produced as b = A*x.
%
% Reference: C. B. Shaw, Jr., "Improvements of the resolution of
% an instrument by numerical solution of an integral equation",
% J. Math. Anal. Appl. 37 (1972), 83-112.

% Check input.
if (n < 1), error('The order n must be positive'), end

% Initialization.
h = pi/(n-1); A = zeros(n,n);

% Compute the matrix A.
cos_ = cos(-pi/2 + [0:(n-1)]*h);
u_ = pi*sin(-pi/2 + [0:(n-1)]*h);
for i=1:n
    for j=i:n
        if (j ~= (n-i+1))
            cos_sum = cos_(i) + cos_(j);
            u = u_(i) + u_(j);
            A(i,j) = (cos_sum*sin(u)/u)^2;
            A(n-j+1,n-i+1) = A(i,j);
        end
    end
    A(i,n-i+1) = (2*cos_(i))^2;
end
A = A + triu(A,1)'; A = A*h;
for i = 1:n
    A(i,1)=A(i,1)/2;
    A(i,n)=A(i,n)/2;
end

% Compute the vectors x and b.
a1 = 2; c1 = 6; t1 =  .8;
a2 = 1; c2 = 2; t2 = -.5;
if (nargout>1)
  x =   a1*exp(-c1*(-pi/2 + [0:(n-1)]'*h - t1).^2) ...
      + a2*exp(-c2*(-pi/2 + [0:(n-1)]'*h - t2).^2);
  b = A*x;
end
