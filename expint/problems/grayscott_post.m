function yprocessed = grayscott_post(y, problem)
% GRAYSCOTT_POST - Postprocessing of numerical solution to the 
%                  Gray-Scott 1D equation. 
%
% SYNOPSIS:
%   yprocessed = grayscott_post(u, problem);
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
%   GRAYSCOTT, EXPGLM

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.4 $  $Date: 2005/10/22 02:50:13 $

% reshape into matrix form, still in fourier space:
uhat = y(1:problem.ND);
vhat = y(problem.ND+1:end);

% Into physical space:
yprocessed = real(ifft(uhat));
