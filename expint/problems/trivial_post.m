function UF = trivial_post(U, problem)
% TRIVIAL_POST -- Trivial or dummy version of a postprocessing function that 
%                 does nothing to the U-value.
%
% SYNOPSIS:
%   UF = trivial_post(U, problem);
%
% PARAMETERS:
%   U        - Numerical solution at some (usually final) time step.
%   problem  - Problem specific parameters. Not used
%
% RETURNS:
%   UF  - The input value `U'.
%
% SEE ALSO:
%   EXPGLM.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.2 $  $Date: 2005/10/22 02:50:13 $

UF = U;
