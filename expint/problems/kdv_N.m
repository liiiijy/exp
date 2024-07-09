function Nr = kdv_N(U, t, problem)
% KDV_N -- Non-linear term of KdV equation.
%
% SYNOPSIS:
%   Nr = kdv_N(U, t, problem);
% 
% DESCRIPTION:
%   Evaluates the non-linear term of the semi-discretised KDV equation
%   at a point `U' in Fourier space.
%
% PARAMETERS:
%   U       - Evaluation point in Fourier space.
%   t       - timepoint.
%   problem - Problem dependent parameters.  Should be defined by KDV.
%
% RETURNS:
%   Nr      - Value of non-linear term at `U'.
%
% SEE ALSO:
%   KDV, EXPGLM

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.2 $  $Date: 2005/10/22 02:50:13 $

Nr = -0.5i * norm(problem.origy0) * problem.k .* fft(real(ifft(U)).^2);

