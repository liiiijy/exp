function f = paratest_LplusN(t, y, options, problem)
% PARATEST_LPLUSN - Standalone PARATEST right hand side to be used with
% ODE15S.
%
% SYNOPSIS:
%   dy = paratest_LplusN(t, y, opts, pr);
% 
% PARAMETERS:
%   t       - Time.
%   y       - State value at time `t'.
%   opts    - ODE solver options
%   problem - Problem specific parameters as defined by PARATEST.
%
% RETURNS:
%   dy   - Derivative/ODE RHS at (t,y).
%
% SEE ALSO:
%   PARATEST, EXPGLM, ODE15S

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.3 $  $Date: 2005/10/22 02:50:13 $

f = full(problem.L * y + feval(problem.N, y, t, problem));
