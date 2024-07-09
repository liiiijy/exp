% TESTGRAYSCOTT2D
%
% Integrates the two dimensional Gray-Scott system.
% Plots data at intermediate and final time points.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.2 $  $Date: 2005/10/22 02:50:14 $

h          = 2;
tspan      = [0, 1e4];
timepoints = 0 : 1e3 : 1e4;
ND         = 128;
alpha      = 0.035;
beta       = 0.060;
length     = 1;
scheme     = 'krogstad';

problem = grayscott2d('ND', ND, 'alpha', alpha, ...
                      'beta', beta, 'length', length);

wantcache('no')
[t, y] = expglm(problem, tspan, h, scheme, timepoints);
[X, Y] = meshgrid(problem.x, problem.y);

for k = 1:numel(t),
   u = real(ifft2(reshape(y(k,        1 : ND^2), ND([1, 1]))));
   v = real(ifft2(reshape(y(k, ND^2 + 1 : end),  ND([1, 1]))));

   figure
   subplot(121), surf(X, Y, u), shading interp, view(2), colorbar
   subplot(122), surf(X, Y, v), shading interp, view(2), colorbar
end
