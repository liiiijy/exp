function N = grayscott3d_N(y, t, problem)
% GRAYSCOTT3D_N - Non-linear term of Gray-Scott equation in 3D.
%
% SYNOPSIS:
%   N = grayscott_N(y, t, problem);
% 
% DESCRIPTION:
%   Evaluates the non-linear term of the semi-discretised
%   Gray-Scott 3d problem.
%
% PARAMETERS:
%   y    - Evaluation point in Fourier space (rolled out to vector)
%   t    - Evaluation point in time. Not used in this function.
%   problem   - Problem dependent parameters.  
%               Should be defined by GRAYSCOTT3D.
%
% RETURNS:
%   N      - Value of non-linear term at `U'.
%
% SEE ALSO:
%   GRAYSCOTT3D

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.3 $  $Date: 2005/10/12 16:28:01 $

ND3 = problem.ND ^ 3;
shape = problem.ND([1, 1, 1]);

% reshape into matrix form, in physical space:
u = ifftn(reshape(y(      1 : ND3), shape));
v = ifftn(reshape(y(ND3 + 1 : end), shape));

Nu = fftn(-u.*v.^2 + problem.alpha * (1 - u));
Nv = fftn( u.*v.^2 - (problem.alpha + problem.beta) .* v);

% reshape into vector form, we are now in fourier space:
N = [Nu(:); Nv(:)];
