function varargout = surfplot(data, timepoints, problem, varargin)
% SURFPLOT
% Plots a solution candidate as a function of time and space. This 
% particular script only works as intended for 1D-problems.
%
% SYNOPSIS
%        surfplot(data, timepoints, problem);
%        surfplot(data, timepoints, problem, 'top');
%
% DESCRIPTION:
%   Plots the data supplied as a function of timepoints supplied 
%   and the x-values found in the problem structure.
%
% PARAMETERS:
%   data - cell-array of solution data. Each cell is the solution at some
%          point in time and is assumed to be a vector (1D only),
%          that is data{k}(:) = u(x, timepoints(k))
%   timepoints - points in time at which the data approximates the
%          some solution.
%   problem - a problem structure. Used for pulling out problemname and
%          x-values.
%   'top'  - the string 'top' if a top view is wanted.
%
% EXAMPLE:a
%  tspan = [0 1]; h = dt(end); h = 0.1; timepoints = [tspan(1):h:tspan(end)];
%  [t, U, UF] = expglm(problem, tspan, h, 'lawson4', timepoints);
%  surfplot(UF, t, problem);

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.4 $  $Date: 2005/10/22 02:50:13 $

nargchk(3,4, nargin);

if min(size(data{1})) ~= 1,
   error('surfplot:only1dproblems', ['This script only supports' ...
       ' 1d-problems']);
end

figure; % Don't overwrite existing figures.
[X, T] = meshgrid(problem.x, timepoints);



% Masssage the incoming datas:
data = cell2mat(data.');
% now data(:,k) is u(x, timespoints(k));

surf(X, T, data, ...
    'FaceColor','interp',...
    'EdgeColor','none',...
    'FaceLighting','phong')
axis tight

if nargin > 3 && strcmp(varargin{1}, 'top'),
   % Top view:
   view(0, 90); 
   colormap(autumn);
   alpha(1);
else
   view(-50,30)
   camlight left
   alpha(0.6);
end

xlabel('x');
ylabel('time');
zlabel('U(x,t)');

title(problem.problemname);
