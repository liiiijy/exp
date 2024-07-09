function N = paratest_N(U, t, problem)
% PARATEST_N -- Non-linear term of PARATEST equation.
%
% SYNOPSIS:
%   Nr = paratest_N(U, problem);
% 
% PARAMETERS:
%   U       - Evaluation point in physical space.
%   problem - Problem dependent parameters.  Defined by PARATEST.
%
% RETURNS:
%   Nr      - Value of non-linear term at `U'.
%
% SEE ALSO:
%   PARATEST, EXPGLM.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.2 $  $Date: 2005/10/22 02:50:13 $

  
dx = problem.spatiallength/problem.ND;

e = ones(problem.ND-1, 1);
z = zeros(problem.ND-1, 1);
N1 = spdiags([-e z e], -1:1, problem.ND-1, problem.ND-1);
D = N1 ./ (2*dx);
DU = D*U;
N = -U .* DU;

x = problem.x;
Uex = x .* (1-x) .* exp(-t);
Uex_x = (1-2*x) .* exp(-t);

PHI = -Uex + 2*exp(-t) + Uex .* Uex_x;

N = N + PHI;
