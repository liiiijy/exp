function uex = schrodingertype_exact(x, t, problem)
% SCHRODINGERTYPE_EXACT - exact solution of the schrodingertype problem
%
% SYNOPSIS:
%   Nr = schrodingertype_exact(x, t, problem);
% 
% PARAMETERS:
%   x       - spatial point
%   t       - point in time
%   problem - Problem dependent parameters.  Defined by SCHRODINGERTYPE.
%
% RETURNS:
%   uex     - the exact value at u(x,t).
%
% SEE ALSO:
%   SCHRODINGERTYPE, EXPGLM.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.2 $  $Date: 2005/10/22 02:50:13 $

uex = x.* (1-x) .* exp(-t);
