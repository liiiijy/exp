function UF = kursiv_post(U, problem)
% KURSIV_POST -- Postprocessing of numerical solution to KURSIV problem
%
% SYNOPSIS:
%   UF = kursiv_post(U, problem);
%
% PARAMETERS:
%   U        - Numerical solution at some (usually final) time step.
%   problem  - Problem specific parameters of the KURSIV problem.  Defined by
%              KURSIV.
%
% RETURNS:
%   UF  - Processed numerical solution `U'.
%
% SEE ALSO:
%   KURSIV, EXPGLM.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.4 $  $Date: 2005/10/22 02:50:13 $

UF = real(ifft(U));
