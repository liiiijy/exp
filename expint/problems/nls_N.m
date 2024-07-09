function Nr = nls_N(U, t, problem)
% NLS_N - Non-linear term of NLS equation.
%
% SYNOPSIS:
%   Nr = nls_N(U, t, problem);
% 
% DESCRIPTION:
%   Evaluates the non-linear term of the semi-discretised non-linear
%   Schrödinger equation at a point `U' in Fourier space.
%
% PARAMETERS:
%   U       - Evaluation point in Fourier space.
%   t       - Evaluation point in time. Not used in this function.
%   problem - Problem dependent parameters.  Should be defined by NLS.
%
% RETURNS:
%   Nr      - Value of non-linear term at `U'.
%
% SEE ALSO:
%   NLS

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.4 $  $Date: 2005/10/13 13:37:40 $

psi = ifft(U);
Nr  = -1i * fft((problem.v + (problem.lambda * psi .* conj(psi))) .* psi);
