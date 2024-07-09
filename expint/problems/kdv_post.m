function UF = kdv_post(U, problem)
% KDV_POST -- Postprocessing of numerical solution to KDV problem
%
% SYNOPSIS:
%   uf = kdv_POST(u, pr);
%
% PARAMETERS:
%   u   - Numerical solution at some (usually final) time step.
%   pr  - Problem specific parameters of the NLS problem.  Defined by
%         KDV.
%
% RETURNS:
%   uf  - Processed numerical solution `u'.
%
% SEE ALSO:
%   KDV, GLOBALORDER

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.2 $  $Date: 2005/04/19 08:33:12 $

UF = real(ifft(U));
