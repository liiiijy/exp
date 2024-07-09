function N = hochost_N(U, t, problem)
% HOCHOST_N -- Non-linear term of problem of Hocbbruck and Ostermann.
%
% SYNOPSIS:
%   Nr = hochost_N(U, problem);
% 
% PARAMETERS:
%   U       - Evaluation point in physical space.
%   problem - Problem dependent parameters.  Defined by HOCHOST.
%
% RETURNS:
%   Nr      - Value of non-linear term at `U'.
%
% SEE ALSO:
%   HOCHOST, EXPGLM.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.3 $  $Date: 2005/10/22 02:50:13 $

N = 1 ./ (1 + U.^2);

x = problem.x;
Uex = x .* (1-x) .* exp(t);

PHI = Uex + 2*exp(t) - 1./(1+Uex.^2);

N = N + PHI;
