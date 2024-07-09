function varargout = orderplot(M, thexlabel, theylabel, thetitle)
% ORDERPLOT
% Produces local or global error plot. Plots local or global error versus
% timestep. Error-data and timestep are read from the M-structure, and
% typically stem from either localorder or globalorder.
%
% SYNOPSIS: 
%         orderplot(M, thexlabel, theylabel, thetitle)
%     h = orderplot(M, thexlabel, theylabel, thetitle)
% 
%
% PARAMETERS:
%     M - cellarray of structures, one for each scheme to be plotted.
%         Typically initially set up with SETUPSCHEMES, then 
%         data is added by LOCALORDER or GLOBALORDER.
%     thexlabel - string, labelling of timestep-axis.
%     theylabel - string, labelling of error-axis.
%     thetitle  - string, title of plot. Typical choice is problem.problemname.
%
% SEE ALSO:
%   LOCALORDER, GLOBALORDER, SETUPSCHEMES, TIMINGPLOT

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.4 $  $Date: 2005/04/28 09:10:59 $

error(nargoutchk(0, 1, nargout));

LineWidth = 2;
Nschemes = length(M);

font = ''; % Used for labels.

h = figure;
axes('FontSize', 14);


figure(h);
hold on;

for k = 1:Nschemes,
   plot(M{k}.hs, M{k}.error, M{k}.linestyle, 'LineWidth', 2);
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
