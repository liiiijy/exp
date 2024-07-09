function schemes = localorder(problem, dt, schemes, varargin)
% LOCALORDER    
%   Does a local order test of exponential integrators
%
% SYNOPSIS:
%    schemes         = localorder(problem, dt, schemes);
%    schemes         = localorder(problem, dt, schemes, exactsolver);
%
% PARAMETERS:
%   problem -
%         A structure describing the problem in question. Contains 
%         the linear part, the non-linear part, initial condition and
%         other relevant information.
%   dt  - Vector of time steps.
%   schemes- cell of structures defining steppers investigated 
%         in experiment. Use SETUPSCHEMES to define this structure and
%         and for documentation.
%   exactsolver - a string, either 'ode15s', or some other expint-scheme
%         available in PATH. If not ode15s, the smallest timestep will 
%         be divided by dtMult=10 for the reference solution computation.
%         If this string is 'exact', then the field exactsolution in the
%         problem will be used.
%
%
% RETURNS:
%   schemes - An augmented version of the incoming 'schemes'-structures, 
%             now with local error results and corresponding timesteps
%             for all schemes. Timing results is also provided, but it
%             is up to the user to interpret them sensibly. Usually, 
%             a global order test is more appropriate for timing.
%
%
% SEE ALSO:
%   GLOBALORDER, SETUPSCHEMES, ORDERPLOT.
%

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.6 $  $Date: 2005/10/22 02:50:13 $

error(nargchk(3, 4, nargin));

% Number of steps pr. min stepsize in refsol.
dtMult = 10; 

% Relevant if ode15s is used:
ode15soptions = odeset('RelTol', 1.0e-13, 'AbsTol', 1.0e-15);

% Check if user has some preference on the exactsolver:
if (nargin > 3),
   exactsolver = varargin{1};
   if (ischar(exactsolver)),
      % Here we only check if the supplied argument is valid, that
      % is, either a function in path, or 'ode15s':
      if ~(exist(exactsolver)==2 || exist(exactsolver)==3 || strcmpi(exactsolver,'ode15s'))
	 error('localorder:noexactsolver', ['Scheme ', exactsolver, ...
		' not found in path']);
      end
   end
else
   % Default exactsolver:
   exactsolver = 'ode15s';
end


disp(['Local order: ', problem.problemname]);


inputexps = numel(dt);
dt = unique(dt);
if (length(dt) ~= inputexps)
   warning('localorder:nonuniquedt', ...
       'There were identical elements in the supplied dt-vector');
end


%% Exact solution. May be supplied in each problem, or else, it
%% will be computed numerically.
[tspan, sortidx] = sort(dt);
if strcmp(exactsolver, 'exact'),
   if isfield(problem, 'exactsolution'),
      disp([' Exact solution provided by problem']);
      refsols = zeros(length(dt) + 1, length(problem.y0));
      refsols(1,:)  = problem.y0.'; % Not used, but added for consistency
      for e=1:length(dt),
	 refsols(e+1, :) = feval(problem.exactsolution, problem.x, dt(e), problem).';
      end
   else
      error('localorder:noexactsolver', ...
	  'Problem does not have an exact solution available');
   end
else
   disp([datestr(clock, 31), ' |-> Starting reference solution...']);
   

   if strcmp(exactsolver, 'ode15s'),
      disp([' Exact solver: ', exactsolver]);      
      t = cputime;
      [ts, refsols] = ode15s(problem.LplusN, [0 tspan], ...
	  problem.y0, ode15soptions, problem);
      t = cputime - t;
      
   else
      % Use an expintscheme for the exact solution, but 
      % with a timestep much smaller than the smallest timestep
      % used in the global order test.
      hmin = min(dt)/dtMult;
      disp([' Exact solver: ', exactsolver, ', h = ', num2str(hmin)]);
      t = cputime;
      [ts, refsols] = expglm(problem, [0 max(tspan)], ...
	  hmin, exactsolver, [0 tspan]);
      t = cputime - t;
   end
  
   disp([datestr(clock, 31), ' |-> Reference solution computed in ', ...
	  num2str(t, '%0.4g'), ' seconds']);
end

% Some massaging of the schemes-cell-structure:
for k = 1:length(schemes),

   % Make space for the error-results:
   schemes{k}.error = zeros(length(dt), 1);
   
   % Make space for timing-results, 
   schemes{k}.timing = zeros(length(dt), 1);
   
   % We insert the timesteps used in the scheme structure as well.
   schemes{k}.hs = dt;

   % Grab out relstages if present for the schemes. 
   % This is if some schemes need smaller timestep in order
   % to be compared fairly.
   if isfield(schemes{k}, 'relstages'),
      relstages(k) = schemes{k}.relstages;
   else
      relstages(k) = 1;
   end
end

cput = cputime;

%% Do experiments for each timestep and each solver:
for e = 1:length(dt), % for each dt
   disp([datestr(clock, 31), ' |-> Experiment ', num2str(e, '%03d'), ...
	  ', dt = ', num2str(dt(e), '%0.5e')]);
   
   for k = 1:length(schemes), % for each scheme
      tic
      [t, numsol] = expglm(problem, [0 dt(e)], ...
	  dt(e) / relstages(k), schemes{k}.coffcn);
      schemes{k}.error(e) = norm(refsols(1+sortidx(e),:) - numsol);
      % We add 1 to the refsols-index because timepoint 0 is included
      % in refsols (= initial condition).
      schemes{k}.timing(e) = toc;
   end
end

cput = cputime - cput;
disp([datestr(clock, 31), ' |-> Experiments completed in ', ...
       num2str(cput, '%0.4g'), ' seconds']);
