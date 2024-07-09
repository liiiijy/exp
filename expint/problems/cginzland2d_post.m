function uf = cginzland2d_post(u, problem)
% CGINZLAND2D_POST
% Postprocessing of numerical solution to the Complex Ginzburg-Landau 2D
% equation. 
%
% SYNOPSIS:
%   uf = cginzland2d_post(u, problem);
%
% PARAMETERS:
%   u   - Numerical solution at some point in time, row-vector of
%         dimension ND^2 (rolled out).
%   problem  - Problem specific parameters.
%
% RETURNS:
%   uf  - Processed numerical solution `u', in matrix form.
%
% SEE ALSO:
%   CGINZLAND2D, EXPGLM

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.2 $  $Date: 2005/10/22 02:50:13 $

% reshape into matrix form, still in fourier space:
uhat = reshape(u, [problem.ND, problem.ND]);

% Into physical space:
uf = ifft2(uhat);
