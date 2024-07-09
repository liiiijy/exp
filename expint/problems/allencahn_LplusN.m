function f = allencahn_LplusN(t, y, options, problem)
% AC_LPLUSN - Standalone Allen--Cahn right hand side to be used with ODE15S.
%
% SYNOPSIS:
%   f = allencahn_LplusN(t, y, opts, pr);
% 
% PARAMETERS:
%   t    - Time.
%   y    - State value at time `t'.
%   opts - ODE solver options
%   pr   - Problem specific parameters as defined by AC.
%
% RETURNS:
%   f   - Derivative/ODE RHS at (t,y).
%
% SEE ALSO:
%   ALLENCAHN, EXPGLM, ODE15S

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.3 $  $Date: 2005/10/22 02:50:13 $

f = problem.L * y + feval(problem.N, y, t, problem);
