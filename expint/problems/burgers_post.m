function UF = burgers_post(U, problem)
% BURGERS_POST -- Postprocessing of numerical solution to BURGERS problem
%
% SYNOPSIS:
%   UF = burgers_POST(U, problem);
%
% PARAMETERS:
%   U        - Numerical solution at some (usually final) time step.
%   problem  - Problem specific parameters of the BURGERS problem.  Defined by
%              BURGERS.
%
% RETURNS:
%   UF  - Processed numerical solution `U'.
%
% SEE ALSO:
%   BURGERS, EXPGLM.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.5 $  $Date: 2005/10/22 02:50:13 $

UF = real(ifft(U));
