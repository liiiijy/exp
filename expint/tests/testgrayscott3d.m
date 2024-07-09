% TESTGRAYSCOTT3D
%
% Integrates the three dimensional Gray-Scott system.
% Plots data at endpoint.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.2 $  $Date: 2005/10/22 02:50:14 $

h      = 1;
tspan  = [0, 4e3];
ND     = 48;
alpha  = 0.025;
beta   = 0.060;
length = 0.75;
scheme = 'modgenlawson45';

problem = grayscott3d('ND', ND, 'alpha', alpha, ...
                      'beta', beta, 'length', length);

wantcache('no')
[t, y] = expglm(problem, tspan, h, scheme);

u = real(ifftn(reshape(y(       1 : ND^3), ND([1, 1, 1]))));
v = real(ifftn(reshape(y(ND^3 + 1 : end),  ND([1, 1, 1]))));
[X, Y, Z] = meshgrid(problem.x, problem.y, problem.z);

figure
uval = 0.3;
ip(1) = patch(isosurface(X, Y, Z, u, uval));
isonormals(X, Y, Z, u, ip(1))
cp(1) = patch(isocaps(X, Y, Z, u, uval), ...
              'FaceColor', 'interp', 'EdgeColor', 'none');
colormap(flipud(hsv(256)))
set(gcf, 'Renderer', 'zbuffer')
camlight, camlight left, lighting phong
view(3), axis tight, camzoom(1.15), colorbar

figure
vval = 0.25;
ip(2) = patch(isosurface(X, Y, Z, v, vval), ...
              'FaceColor', 'blue', 'EdgeColor', 'none');
isonormals(X, Y, Z, v, ip(2))
cp(2) = patch(isocaps(X, Y, Z, v, vval), ...
              'FaceColor', 'interp', 'EdgeColor', 'none');
colormap(flipud(hsv(256)))
set(gcf, 'Renderer', 'zbuffer')
camlight, camlight left, lighting phong
view(3), axis tight, camzoom(1.15), colorbar
