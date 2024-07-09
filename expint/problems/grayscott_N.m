function N = grayscott_N(y, t, problem)
% GRAYSCOTT_N - Non-linear term of Gray-Scott equation in 1D.
%
% SYNOPSIS:
%   N = grayscott_N(y, t, problem);
% 
% DESCRIPTION:
%   Evaluates the non-linear term of the semi-discretised
%   Gray-Scott 1d problem.
%
% PARAMETERS:
%   y       - Evaluation point in Fourier space (rolled out to vector)
%   t       - Evaluation point in time. Not used in this function.
%   problem - Problem dependent parameters. Should be defined by GRAYSCOTT.
%
% RETURNS:
%   N      - Value of non-linear term at `U'.
%
% SEE ALSO:
%   GRAYSCOTT

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.2 $  $Date: 2005/10/12 16:28:01 $

% reshape into matrix form, in physical space:
u = ifft(y(1 : problem.ND));
v = ifft(y(problem.ND + 1 : end));

Nu = fft(-u.*v.^2 + problem.alpha * (1 - u));
Nv = fft( u.*v.^2 - (problem.alpha + problem.beta) .* v);

% reshape into vector form, we are now in fourier space:
N = [Nu(:); Nv(:)];
