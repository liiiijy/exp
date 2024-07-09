function N = hochostint_N(U, t, p)
% HOCHOSTINT_N -- Non-linear term of parabolic equation used by
% Hochbruck and Ostermann, integral version.
%
% SYNOPSIS:
%   Nr = hochostint_N(U, t, problem);
%
% PARAMETERS:
%   U       - Evaluation point in physical space.
%   t       - Time.
%   problem - Problem dependent parameters.  Defined by HOCHOSTINT.
%
% RETURNS:
%   Nr      - Value of non-linear term at `(t, U)'.
%
% SEE ALSO:
%   HOCHOSTINT, EXPGLM

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.4 $  $Date: 2005/10/22 02:50:13 $

x = p.x;
dx = p.spatiallength / p.ND;

% Approximate value of \int_0^1 u(x,t) dx 
% Note:
%   This implicitly relies on homogeneous Dirichlet boundary conditions
%   on both boundaries.  Changing boundary conditions requires changing
%   the call to SUM.

if mod(p.ND, 2) == 0,
   % Simpson approx
   s = 4*sum(U(1:2:end)) + 2*sum(U(2:2:end));
   s = s * dx / 3;
else
   % Trapezoidal rule
   s = dx * sum(U);
end

N = (x .* (1 - x) + 11/6).*exp(t) + s;
