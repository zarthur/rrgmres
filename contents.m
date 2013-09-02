% Range Restricted GMRES Demo.
% Version 1.0  May 2011.
% Copyright (c) 2011.
%
% Demonstration.
%   rrgmres_demo.fig        -   The graphical interface for the RRGMRES 
%                               demo for linear discrete ill-posed 
%                               problems with a square nonsymmetric 
%                               matrix.
%
%   rrgmres_demo.m          -   The RRGMRES demo for linear discrete 
%                               ill-posed problem with a square 
%                               nonsymmetric matrix.
%
%   sym_rrgmres_demo.fig    -   The graphical interface for demo for the
%                               RRGMRES method for linear discrete 
%                               ill-posed problems with a symmetric 
%                               matrix.
%
%   sym_rrgmres_demo.m      -   The demo for the RRGMRES method for linear 
%                               discrete ill-posed problems with a 
%                               symmetric matrix.
%
% Algorithms
%   rrgmres_dp.m            -   The RRGMRES algorithm for linear discrete 
%                               ill-posed problems with a square 
%                               nonsymmetric matrix. This version uses the 
%                               discrepancy principle to decide when to 
%                               terminate the iterations.
%
%   rrgmres_iter.m          -   The RRGMRES algorithm for linear discrete 
%                               ill-posed problems with a square
%                               nonsymmetric Matrix. This version allows 
%                               the user to specify the desired number of 
%                               iterations.
%
%   sym_rrgmres_dp.m        -   The RRGMRES algorithm for linear discrete 
%                               ill-posed problems with a symmetric 
%                               matrix. This version uses the discrepancy 
%                               principle to decide when to terminate the 
%                               iterations.
%
%   sym_rrgmres_iter.m      -   The RRGMRES algorithm for linear discrete 
%                               ill-posed problems with a symmetric 
%                               matrix. This version allows the user to 
%                               specify the number of desired iterations.
%
% Test problems.
%   baart_alt.m             -   Discretization of the Fredholm integral 
%                               equation of the first kind described by 
%                               Baart using a Nystrom method based on the 
%                               composite trapezoidal rule with 
%                               equidistant nodes. The linear discrete 
%                               ill-posed problem obtained has a square 
%                               nonsymmetric matrix.
%
%   deriv2_alt.m            -   Discretization of a Fredholm integral 
%                               equation of the first kind that is a 
%                               Green’s function for the second derivative 
%                               on the interval [0, 1]. The discrete 
%                               problem is obtained by a Nystrom method 
%                               based on the composite trapezoidal rule 
%                               with a square matrix. The linear discrete 
%                               ill-posed problem obtained has a square 
%                               nonsymmetric matrix. The solution can be 
%                               chosen to be a discretized linear or 
%                               exponential function.
%
%   phillips_alt.m          -   Discretization of the Fredholm integral 
%                               equation of the first kind described by 
%                               Phillips using a Nystrom method based on 
%                               the composite trapezoidal rule with 
%                               equidistant nodes. The linear discrete 
%                               ill-posed problem obtained has a square 
%                               nonsymmetric matrix.
%
%   shaw_alt.m              -   Discretization of the Fredholm integral of 
%                               the first kind discussed by Shaw [9] using 
%                               a Nystrom method based on the composite 
%                               trapezoidal rule with equidistant nodes. 
%                               The linear discrete ill-posed problem 
%                               obtained has a square nonsymmetric matrix.
%
% Further information is provided in primer.pdf
