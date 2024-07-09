function uex = hochost_exact(x, t, problem)
% HOCHOST_EXACT - exact solution of the hochost problem
%
% SYNOPSIS:
%   Nr = hochost_exact(x, t, problem);
% 
% PARAMETERS:
%   x       - spatial point
%   t       - point in time
%   problem - Problem dependent parameters.  Defined by HOCHOST.
%
% RETURNS:
%   uex     - the exact value at u(x,t).
%
% SEE ALSO:
%   HOCHOST, EXPGLM.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.3 $  $Date: 2005/10/22 02:50:13 $

uex = x.* (1-x) .* exp(t);
