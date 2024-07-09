%
% A sample test file to test global order behaviour for some exponential
% integrators
%

% This file is part of the 'Expint'-package, 
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.1 $ $Date: 2005/04/19 08:45:35 $

% Set up choice of schemes and placeholder for results:
M = setupschemes('abnorsett4', 'lawson4', 'hochost4', ...
                 'genlawson43', 'rkmk4t', 'cfree4', ...
                 'cranknicolson');

% Set relative stages for the first method.
M{1}.relstages = 4;

% Set the problems and user defined data.
problem1 = kursiv('ND', 128, 'IC', 'smooth');
tspan1 = [0, 65];
steplist1 = round(10.^(6.5:-0.05:2)); 
dt1 = unique(diff(tspan1)./steplist1);

problem2 = kdv('ND', 128, 'IC', 'sechsqr');
tspan2 = [0, 2*pi/625];
steplist2 = round(10.^(6.5:-0.05:0)); 
dt2 = unique(diff(tspan2)./steplist2);

problem3 = ac('ND', 50, 'lambda', 0.001, 'IC', 'smooth');
tspan3 = [0, 3];
steplist3 = round(10.^(5:-0.05:0)); 
dt3 = unique(diff(tspan3)./steplist3);

problem4 = hoInt('ND', 100, 'IC', 'smooth');
tspan4 = [0, 1];
steplist4 = round(10.^(4:-0.02:0)); 
dt4 = unique(diff(tspan4)./steplist4);

problem5 = nls('ND', 512, 'lambda', 1, 'IC', 'smooth', 'Pot', 'reg2');
tspan5 = [0, 1];
steplist5 = round(10.^(4:-0.02:0)); 
dt5 = unique(diff(tspan5)./steplist5);

% Do the global order test. The M structure is returned augmented.
M1 = globalorder(problem1, tspan1, dt1, M, 'hochost4');
M2 = globalorder(problem2, tspan2, dt2, M, 'hochost4');
M3 = globalorder(problem3, tspan3, dt3, M, 'abnorsett4');
M4 = globalorder(problem4, tspan4, dt4, M, 'abnorsett4');
M5 = globalorder(problem5, tspan5, dt5, M, 'hochost4');


% Plot the results
orderplot(M1, 'Timestep h', 'Error', ['Global error, ', ...
                    problem1.problemname]);
orderplot(M2, 'Timestep h', 'Error', ['Global error, ', ...
                    problem2.problemname]);
orderplot(M3, 'Timestep h', 'Error', ['Global error, ', ...
                    problem3.problemname]);
orderplot(M4, 'Timestep h', 'Error', ['Global error, ', ...
                    problem4.problemname]);
orderplot(M5, 'Timestep h', 'Error', ['Global error, ', ...
                    problem5.problemname]);



