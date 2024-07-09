function Nr = schrodingertype_N(U, t, problem)
% SCHRODINGERTYPE_N - Non-linear term of SCHRODINGERTYPE equation.
%
% SYNOPSIS:
%   Nr = schrodingertype_N(U, t, problem);
% 
% DESCRIPTION:
%   Evaluates the non-linear term of the semi-discretised non-linear
%   Schrödinger equation at a point `U' in Fourier space.
%
% PARAMETERS:
%   U       - Evaluation point in Fourier space.
%   t       - Evaluation point in time. Not used in this function.
%   problem - Problem dependent parameters.  Should be defined by SCHRODINGERTYPE.
%
% RETURNS:
%   Nr      - Value of non-linear term at `U'.
%
% SEE ALSO:
%   SCHRODINGERTYPE

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.2 $  $Date: 2005/10/13 13:37:40 $

dx = problem.spatiallength/problem.ND;
 
e = ones(problem.ND-1, 1);
z = zeros(problem.ND-1, 1);
N1 = spdiags([-e z e], -1:1, problem.ND-1, problem.ND-1);
D = N1 ./ (2*dx);
DU = D*U;
N  = 1i * (U .* DU);

x = problem.x;
Uex  = x .* (1-x) .* exp(-t);
Uex_x = (1-2*x) .* exp(-t);

iPHI = Uex + 2i*exp(-t) + 1i * Uex .* Uex_x; 

Nr = N - iPHI;
