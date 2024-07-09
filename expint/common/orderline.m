function orderline(gradient, varargin)
% ORDERLINE
% Plots a straight line in the current plot. Useful for visualising
% which order a specific numerical schemes has in an order-plot.
% May be used both interactive and non-interacive.
% 
%
% SYNOPSIS:
%   orderline(gradient)
%   orderline(gradient, [h, err])
%   orderline(gradient, [h, err], textplacement)
%
% PARAMETERS:
%   gradient - the derivative of the straigt line (the order)
%   [h, err] - one point in the plot which the line is to go through
%              if omitted, ginput is used.
%   textplacement - a factor used for placing the gradient as text
%                   in the plot, DEFAULT 0.28.
%

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.2 $ $Date: 2005/04/18 11:57:15 $

% To be edited by the user of the script:
hs = [10^-5 10^-0.01];
% This is the time-step interval for the order-line. This
% is best to edit manually instead of trying to calculate
% automatically.

% Save current axis settings:
oldaxes = axis; 

% Interactive use if only one input argument:
if (nargin == 1) 
  [h, err] = ginput(1); % Gets coordinates from mouse.
  disp(['[h, err] = [', num2str(h), ', ', num2str(err), ']']);
end

if (nargin > 1)
  point = varargin{1};
  h = point(1);
  err = point(2);
end

if (nargin > 2)
  textplacement = varargin{2};
else
  textplacement = 0.28;
end

% The b in "y = ax+b":
b = log10(err) - gradient*log10(h);

% Draw the line:
loglog(hs, 10.^(gradient*log10(hs)+b), 'k:', 'LineWidth', 2)

% Place the steepness/order/gradient in the plot as text:
text(h,textplacement*err, sprintf('%g', gradient))
  
axis(oldaxes);
