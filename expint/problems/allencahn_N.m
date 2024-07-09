function Nr = allencahn_N(u, t, problem)
% ALLENCAHN_N -- Non-linear term of Allen-Cahn equation.
%
% SYNOPSIS:
%   Nr = allencahn_N(u, t, problem);
% 
% DESCRIPTION:
%   Evaluates the non-linear term of the semi-discretised Allen--Cahn
%   equation at a point `u' in physical space. To deal with boundary
%   conditions, we work with y = x + u, where u is dealt with with
%   homogeneous BCs.
%
% PARAMETERS:
%   u       - Evaluation point in physical space.
%   t       - Time. 
%   problem - Problem dependent parameters.  Defined by ALLENCAHN.
%
% RETURNS:
%   Nr      - Value of non-linear term at `u'.
%
% SEE ALSO:
%   ALLENCAHN, EXPGLM.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.4 $  $Date: 2005/10/22 02:50:13 $

Nr = (problem.x + u) - (problem.x + u).^3;
