%
% A sample test file to test global order behaviour for some exponential
% integrators
%

% This file is part of the 'Expint'-package, 
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.1 $ $Date: 2005/05/11 08:56:49 $

problem = kursiv('ND', 32);

timespan = [0, 1];
steplist = round(10.^(3:-0.1:0)); 
dt = diff(timespan)./steplist;

% Set up choice of schemes and placeholder for results:
M = setupschemes('lawson4', 'hochost4', 'etd4rk', 'rkmk4t',  ...
    'abnorsett4', 'ablawson4', 'etd5rkf', ...
    'genlawson45', 'modgenlawson45' );


% Do the global order test. The M structure is returned augmented.
M = globalorder(problem, timespan, dt, M);

% Plot the results
orderplot(M, 'Timestep h', 'Global error', problem.problemname);
timingplot(M, 'Time used', 'Global error', problem.problemname);
