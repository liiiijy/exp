function yprocessed = grayscott3d_post(y, problem)
% GRAYSCOTT3D_POST - Postprocessing of numerical solution to the 
%                    Gray--Scott 3D equation. 
%
% SYNOPSIS:
%   yprocessed = grayscott3d_post(y, problem);
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
%   GRAYSCOTT3D, EXPGLM

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.3 $  $Date: 2005/10/22 02:50:13 $

% reshape into matrix form, still in fourier space:
uhat = reshape(y(1:problem.ND^3), [problem.ND, problem.ND, problem.ND]);
vhat = reshape(y(problem.ND^3+1:end), [problem.ND, problem.ND, problem.ND]);

% Into physical space:
yprocessed = real(ifftn(uhat));
