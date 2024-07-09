function f = hochost_LplusN(t, y, options, problem)
% HOCHOST_LPLUSN - Standalone HOCHOST right hand side to be used with
%                  ODE15S.
%
% SYNOPSIS:
%   dy = hochost_LplusN(t, y, opts, pr);
% 
% PARAMETERS:
%   t       - Time.
%   y       - State value at time `t'.
%   opts    - ODE solver options
%   problem - Problem specific parameters as defined by HOCHOST.
%
% RETURNS:
%   dy   - Derivative/ODE RHS at (t,y).
%
% SEE ALSO:
%   HOCHOST, EXPGLM, ODE15S

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.5 $  $Date: 2005/10/22 02:50:13 $

f = full(problem.L * y + feval(problem.N, y, t, problem));
