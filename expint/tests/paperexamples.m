%
% The file used to plot the example figures in this package's 
% accompanying paper.
%

% This file is part of the 'Expint'-package, 
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.3 $ $Date: 2005/10/22 02:50:13 $


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Global error and timing plot.

problem = kursiv('ND', 256);

timespan = [0, 1];
steplist = unique(round(10.^(3:-0.01:1)));
dt = diff(timespan)./steplist;

% Set up choice of schemes and placeholder for results:
M = setupschemes('lawson4', 'hochost4', 'etd4rk', 'rkmk4t',  ...
    'abnorsett4', 'ablawson4', 'etd5rkf', ...
    'genlawson45', 'modgenlawson45' );


% Do the global order test. The M structure is returned augmented.
M = globalorder(problem, timespan, dt, M);

% Plot the results
orderplot(M, 'Timestep h', 'Global error', problem.problemname);

print -depsc2 'globalerror.eps';


timingplot(M, 'Time used', 'Global error', problem.problemname);
print -depsc2 'timing.eps';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The example surface plot.

problem = kursiv('ND', 256);
[t, U, UF] = expglm(problem, [0 100], 0.1, 'lawson4', [0:0.2:100]);
surfplot(UF, t, problem);
print -dpng -r300 'foo.png'
