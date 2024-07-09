function f = diagLplusN(t, y, options, problem)
% DIAGLPLUSN - Right hand side for problems suited for Matlabs ODE suite.
% This LplusN is suited for all diagonal problems at least.
%
% SYNOPSIS:
%   dy = diagLplusN(t, y, opts, pr);
% 
% PARAMETERS:
%   t    - Time.
%   y    - State value at time `t'.
%   opts - ODE solver options
%   pr   - Problem structure
%
% RETURNS:
%   dy   - Derivative/ODE RHS at (t,y).
%
% SEE ALSO:
%   ODE15S, EXPGLM

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.3 $  $Date: 2005/10/22 02:50:13 $

f = problem.L .* y + feval(problem.N, y, t, problem);
