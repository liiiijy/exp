function N = grayscott2d_N(y, t, problem)
% GRAYSCOTT2D_N - Non-linear term of Gray-Scott equation in 2D.
%
% SYNOPSIS:
%   N = grayscott2d_N(y, t, problem);
% 
% DESCRIPTION:
%   Evaluates the non-linear term of the semi-discretised
%   Gray-Scott 2d problem.
%
% PARAMETERS:
%   y       - Evaluation point in Fourier space (rolled out to vector)
%   t       - Evaluation point in time. Not used in this function.
%   problem - Problem dependent parameters. Should be defined by GRAYSCOTT2D.
%
% RETURNS:
%   N      - Value of non-linear term at `U'.
%
% SEE ALSO:
%   GRAYSCOTT2D

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.6 $  $Date: 2005/10/12 16:28:01 $

ND2 = problem.ND ^ 2;
shape = problem.ND([1, 1]);

% reshape into matrix form, in physical space:
u = ifft2(reshape(y(      1 : ND2), shape));
v = ifft2(reshape(y(ND2 + 1 : end), shape));

Nu = fft2(-u.*v.^2 + problem.alpha * (1 - u));
Nv = fft2( u.*v.^2 - (problem.alpha + problem.beta) .* v);

% reshape into vector form, we are now in fourier space:
N = [Nu(:); Nv(:)];
