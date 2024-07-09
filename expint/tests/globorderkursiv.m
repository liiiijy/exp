% 
% Code used to reproduce Figure 1 in accompanying paper.
%
%

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.2 $  $Date: 2005/11/08 19:46:22 $

problem = kursiv('ND', 256);
timespan = [0, 1];
steplist = unique(round(10.^(3:-0.01:1)));
dt = diff(timespan)./steplist;

% Set up choice of schemes and placeholder for results:
M = setupschemes('lawson4', 'hochost4', 'etd4rk', 'rkmk4t', ...
                 'abnorsett4', 'ablawson4', 'etd5rkf',      ...
                 'genlawson45', 'modgenlawson45' );

% Do the global order test. The M structure is returned augmented.
M = globalorder(problem, timespan, dt, M, 'ode15s');

% Plot the results
orderplot(M, 'Timestep h', 'Global error', problem.problemname);
timingplot(M, 'Time used', 'Global error', problem.problemname);
