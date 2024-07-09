function UF = nls_post(U, problem)
% NLS_POST -- Postprocessing of numerical solution to NLS problem
%
% SYNOPSIS:
%   uf = nls_post(U, problem);
%
% PARAMETERS:
%   U        - Numerical solution at some (usually final) time step.
%   problem  - Problem specific parameters of the NLS problem.  Defined by
%              NLS.
%
% RETURNS:
%   UF  - Processed numerical solution `U'.
%
% SEE ALSO:
%   NLS, EXPGLM.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.5 $  $Date: 2005/10/22 02:50:13 $

% Return the density of the psi-function according to the common
% interpretation of the solution to the Schrodinger-equation.
UF = abs(ifft(U)) .^ 2;
