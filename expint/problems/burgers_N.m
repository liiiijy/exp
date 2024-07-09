function Nr = burgers_N(U, t, problem)
% BURGERS_N -- Non-linear term of the Burgers equation.
%
% SYNOPSIS:
%   Nr = burgers_N(U, problem);
% 
% DESCRIPTION:
%   Evaluates the non-linear term of the semi-discretised
%   Burgers equation at a point `U'.
%
% PARAMETERS:
%   U       - Evaluation point.
%   problem - Problem dependent parameters.  Defined by BURGERS.
%
% RETURNS:
%   Nr      - Value of non-linear term at `U'.
%
% SEE ALSO:
%   BURGERS, EXPGLM.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.5 $  $Date: 2005/10/22 02:50:13 $

Nr = -0.5i * problem.k .* fft(real(ifft(U)) .^ 2);
