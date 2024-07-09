function sinu = sinegordon_N(y, t, problem);
% SINEGORDON_N - Non-linear term of sine-Gordon equation.
%
% SYNOPSIS:
%   Nr = sinegordon_N(y, t, problem);
% 
% DESCRIPTION:
%   Evaluates the non-linear term of the semi-discretised non-linear
%   Schrödinger equation at a point `y' in Fourier space.
%
% PARAMETERS:
%   y       - Evaluation point in Fourier space.
%   t       - Evaluation point in time. Not used in this function.
%   problem - Problem dependent parameters.  Should be defined by SINEGORDON.
%
% RETURNS:
%   Nr      - Value of non-linear term at `y'.
%
% SEE ALSO:
%   SINEGORDON

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.3 $  $Date: 2005/11/09 03:18:29 $

ND = problem.ND;
sinu = zeros(size(y));
sinu(ND+1:2*ND) = fft( - sin(ifft(y(1:ND))));

