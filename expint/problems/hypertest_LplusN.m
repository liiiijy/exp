function f = hypertest_LplusN(t, y, options, problem)
% HYPERTEST_LPLUSN - Standalone HYPERTEST right hand side to be used with
%                  ODE15S.
%
% SYNOPSIS:
%   dy = hypertest_LplusN(t, y, opts, pr);
% 
% PARAMETERS:
%   t       - Time.
%   y       - State value at time `t'.
%   opts    - ODE solver options
%   problem - Problem specific parameters as defined by HYPERTEST.
%
% RETURNS:
%   dy   - Derivative/ODE RHS at (t,y).
%
% SEE ALSO:
%   HYPERTEST, ODE15S

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.1 $  $Date: 2005/10/13 02:08:49 $

f = full(problem.L * y + feval(problem.N, y, t, problem));
