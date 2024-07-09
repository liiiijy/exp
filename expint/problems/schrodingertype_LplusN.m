function f = schrodingertype_LplusN(t, y, options, problem)
% SCHRODINGERTYPE_LPLUSN - Standalone SCHRODINGERTYPE right hand side to be used with
%                  ODE15S.
%
% SYNOPSIS:
%   dy = schrodingertype_LplusN(t, y, opts, pr);
% 
% PARAMETERS:
%   t       - Time.
%   y       - State value at time `t'.
%   opts    - ODE solver options
%   problem - Problem specific parameters as defined by SCHRODINGERTYPE.
%
% RETURNS:
%   dy   - Derivative/ODE RHS at (t,y).
%
% SEE ALSO:
%   SCHRODINGERTYPE, ODE15S

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.1 $  $Date: 2005/06/24 12:20:53 $

f = full(problem.L * y + feval(problem.N, y, t, problem));
