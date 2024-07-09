function yprocessed = grayscott2d_post(y, problem)
% GRAYSCOTT2D_POST - Postprocessing of numerical solution to the 
%                    Gray--Scott 2D equation. 
%
% SYNOPSIS:
%   yprocessed = grayscott2d_post(u, problem);
%
% PARAMETERS:
%   y   - Numerical solution at some point in time, row-vector of
%         dimension ND^2 (rolled out).
%   problem  - Problem specific parameters.
%
% RETURNS:
%   yprocessed  - Processed numerical solution `u', in matrix form.
%
% SEE ALSO:
%   GRAYSCOTT2D, EXPGLM

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.5 $  $Date: 2005/10/22 02:50:13 $

% reshape into matrix form, still in fourier space:
uhat = reshape(y(1:problem.ND^2), [problem.ND, problem.ND]);
vhat = reshape(y(problem.ND^2+1:end), [problem.ND, problem.ND]);

% Into physical space:
yprocessed = real(ifftn(uhat));
