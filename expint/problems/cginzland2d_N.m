function N = cginzland2d_N(y, t, problem)
% CGINZLAND2D_N -- Non-linear term of complex Ginzburg-Landau
% equation in 2D.
%
% SYNOPSIS:
%   N = cginzland2d_N(y, t, problem);
% 
% DESCRIPTION:
%   Evaluates the non-linear term of the semi-discretised
%   complex Ginzburg-Landau 2D problem.
%
% PARAMETERS:
%   y    - Evaluation point in Fourier space (rolled out to vector)
%   t    - Evaluation point in time. Not used in this function.
%   problem - Problem dependent parameters.  
%             Should be defined by CGINZLAND2D.
%
% RETURNS:
%   N      - Value of non-linear term at `U'.
%
% SEE ALSO:
%   CGINZLAND2D

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.3 $  $Date: 2005/10/12 16:28:01 $

% First into physical space:
u = ifft2(reshape(y, problem.ND([1, 1])));

N = - reshape((1 + 1.0i*problem.beta) * fft2(u .* conj(u) .* u), size(y));

