%
% A sample test file to test global order behaviour for some exponential
% integrators
%

% This file is part of the 'Expint'-package, 
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.1 $ $Date: 2005/04/18 08:33:35 $

problem = burgers('ND', 64, 'lambda', 0.03, 'IC', 'smooth');

tspan = [0, 3];
steplist=round(10.^(3:-0.1:1)); 
dt = diff(tspan)./steplist;

% Set up choice of schemes and placeholder for results:
M = setupschemes('abnorsett4', 'lawson4', 'etd4rk', 'hochost4', 'cranknicolson', ...
    'genlawson43');

% Set relative stages for the first method.
M{1}.relstages = 4;

% Do the global order test. The M structure is returned augmented.
M = globalorder(problem, tspan, dt, M, 'lawson4');

% Plot the results
orderplot(M, 'Timestep h', 'Error', ['Global error, ', problem.problemname]);

