%
% TESTCGINZLAND2D
%
% Integrates the Complex Ginzburg-Landau 2D equation
% and makes the data available in a surface plot.
% 
% Animation is simulated through the 'pause'-command,
% hold down some key to "animate".

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.4 $ $Date: 2005/10/22 02:50:14 $


problem = cginzland2d('ND', 32);
h = 0.1;
displayedh = 1;


disp(problem.problemname);

disp(' * Press a key to display next plot, Ctrl-C to stop');
time = 0;
while true, 
   [t, u, uf] = expglm(problem, [time, time + displayedh], ...
       h, 'hochost4');
   time = time + displayedh;
   disp(['Time = ', num2str(time)]);
   
   
   surf(abs(uf{1})); 
   view(0,90);
   shading interp;
   axis tight
   drawnow
   pause

   % Set initial condition for next call to expglm
   problem.y0 = u.';
end
