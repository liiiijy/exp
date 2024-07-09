function varargout = timingplot(M, thexlabel, theylabel, thetitle)
% TIMINGPLOT
% Produces a plot of error vs timeusage for schemes included in the
% M-structure. The globalorder-script appends timing-results to the
% M-structure which can be used by this script.
%
% SYNOPSIS: 
%       timingplot(M, thexlabel, theylabel, thetitle);
%   h = timingplot(M, thexlabel, theylabel, thetitle);
%
% PARAMETERS:
%     M - cellarray of structures, one for each scheme to be plotted.
%         Typically initially set up with SETUPSCHEMES, then 
%         data is added by LOCALORDER or GLOBALORDER.
%    thexlabel - string with the x-axis label.
%    theylabel - string with the y-axis labe.
%    thetitle  - string with the plot-title.
%
% RETURNS:
%            h - figure handle to the produced plot. Optional.
% SEE ALSO:
%     GLOBALORDER, ORDERPLOT, SETUPSCHEMES, ORDERPLOT
%

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.2 $  $Date: 2005/04/22 08:40:40 $

error(nargoutchk(0, 1, nargout));

LineWidth = 2;

Nschemes = length(M);

font = '';     % Used for labels.

h = figure;
axes('FontSize', 14);


figure(h);
hold on;

for k = 1:Nschemes,
   plot(M{k}.timing, M{k}.error, M{k}.linestyle, 'LineWidth', 2);
end


grid on;

for k=1:length(M),
   legendstrings{k} = M{k}.name;
end
legend(legendstrings(:), 4);

set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

% We do not want minor grid lines to be visible.
grid minor

xlabel(sprintf('%s %s', font, thexlabel));
ylabel(sprintf('%s %s', font, theylabel));
title( sprintf('%s %s', font, thetitle));

if nargout > 0,
   varargout{1} = h;
end
