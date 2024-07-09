%
% A sample test file to test local order behaviour for some exponential
% integrators
%

% This file is part of the 'Expint'-package, 
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.1 $ $Date: 2005/05/11 09:19:00 $

problem = nls('ND', 128, 'IC', 'reg2');

steplist = round(10.^(3:-0.1:0)); 
dt = unique(10^-3 .* steplist);

% Set up choice of schemes and placeholder for results:
M = setupschemes('lawson4', 'hochost4', 'etd4rk', 'rkmk4t',  ...
    'etd5rkf', 'friedli', 'strehmelweiner');

% Do the global order test. The M structure is returned augmented.
M = localorder(problem, dt, M, 'hochost4');

% Plot the results
orderplot(M, 'Timestep h', 'Local error', problem.problemname);
