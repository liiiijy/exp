function Nr = hypertest_N(U, t, problem)
% HYPERTEST_N - Non-linear term of HYPERTEST equation.
%
% SYNOPSIS:
%   Nr = hypertest_N(U, t, problem);
% 
% DESCRIPTION:
%   Evaluates the non-linear term of the semi-discretised non-linear
%   Schrödinger equation at a point `U' in Fourier space.
%
% PARAMETERS:
%   U       - Evaluation point in Fourier space.
%   t       - Evaluation point in time. Not used in this function.
%   problem - Problem dependent parameters.  Should be defined by HYPERTEST.
%
% RETURNS:
%   Nr      - Value of non-linear term at `U'.
%
% SEE ALSO:
%   HYPERTEST

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.2 $  $Date: 2005/10/13 13:37:40 $

N = 1i ./ (1 + U.^2);

x = problem.x;
Uex  = x .* (1-x) .* exp(-t);

iPHI = Uex + 2i*exp(-t) + 1i./(1+Uex.^2);

Nr = N - iPHI;
