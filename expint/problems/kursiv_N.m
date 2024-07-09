function Nr = kursiv_N(U, t, problem)
% KURSIV_N -- Non-linear term for the Kuramoto-Sivashinsky equation.
%
% SYNOPSIS:
%   Nr = kursiv_N(U, problem);
% 
% DESCRIPTION:
%   Evaluates the non-linear term of the semi-discretised
%   Kuramoto-Sivashinsky equation at a point `U'.
%
% PARAMETERS:
%   U       - Evaluation point, in Fourier space.
%   problem - Problem dependent parameters.  Defined by KURSIV.
%
% RETURNS:
%   Nr      - Value of non-linear term at `U'.
%
% SEE ALSO:
%   KURSIV.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.2 $  $Date: 2005/10/12 16:28:01 $

Nr = -0.5i * problem.k .* fft(real(ifft(U)) .^ 2);
