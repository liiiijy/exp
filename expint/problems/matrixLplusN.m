function f = matrixLplusN(t, y, options, problem)
% MATRIXLPLUSN - Right hand side for problems suited for Matlabs ODE suite.
% This LplusN is suited for problems where the L-matrix is not
% diagonal, use diagLplusN for diagonal problems.
%
% SYNOPSIS:
%   dy = matrixLplusN(t, y, opts, pr);
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
%   ODE15S, EXPGLM, DIAGLPLUSN

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.2 $  $Date: 2005/10/22 02:50:13 $

f = problem.L * y + feval(problem.N, y, t, problem);
