function schemes = globalorder(problem, tspan, dt, schemes, varargin)
% GLOBALORDER    - Does a global order test of exponential integrators
%
% SYNOPSIS:
%    schemes = globalorder(problem, tspan, dt, schemes);
%    schemes = globalorder(problem, tspan, dt, schemes, exactsolver);
%
% PARAMETERS:
%   problem -
%         A structure describing the problem in question. Contains
%         the linear part, the non-linear part, initial condition and
%         other relevant information.
%   tspan -
%         Time span of integration.  Modelled on similar parameter of
%         ODE45.
%   dt  - Vector of time steps.
%   schemes- cell of structures defining steppers investigated
%         in experiment. Use SETUPSCHEMES to define this structure and
%         and for documentation.
%   exactsolver - a string, either 'ode15s', or some other expint-scheme
%         available in PATH. If not ode15s, the smallest timestep will
%         be divided by dtMult=10 for the reference solution computation.
%         If string is 'exact', the field exactsolution in the problem
%         will be used.
%
% RETURNS:
%   schemes - An augmented version of the incoming 'schemes'-structures,
%             now with extra fields added with data from the experiments
%             usable for plotting or other things.
%
% SEE ALSO:
%   ODE15S, SETUPSCHEMES, ORDERPLOT

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.13 $  $Date: 2005/10/22 02:50:13 $


error(nargchk(4, 5, nargin));
error(nargoutchk(1, 1, nargout));

% Number of steps pr. min stepsize in refsol.
dtMult = 10;

% Relevant if ode15s is used:
ode15soptions = odeset('RelTol', 5.0e-14, 'AbsTol', 1.0e-15);

% Number of experiments
inputexps = numel(dt);
dt = unique(dt); Nexp = numel(dt);
if (length(dt) ~= inputexps)
   warning('globalorder:nonuniquedt', ...
       'There were identical elements in the supplied dt-vector');
end

% Check if user has some preference on the exactsolver:
if nargin > 4,
   exactsolver = varargin{1};

   if ischar(exactsolver),
      % Here we only check if the supplied argument is valid, that
      % is, either a function in path, or 'ode15s' or 'exact':
      if ~(exist(exactsolver)==2 || exist(exactsolver)==3 || ...
          strcmpi(exactsolver, 'ode15s') || strcmpi(exactsolver, 'exact')),
         error('globalorder:noexactsolver', ...
                ['Scheme ', exactsolver, ' not in path']);
      end
   end
else
   % Default exactsolver:
   exactsolver = 'ode15s';
end

disp(['Global order: ', problem.problemname]);


%% Exact solution. May be supplied in each problem, or else, it
%% will be computed numerically.


if strcmp(exactsolver, 'exact'),
   if isfield(problem, 'exactsolution'),
      disp([' Exact solution provided by problem']);
      refsol = feval(problem.exactsolution, problem.x, tspan(end), problem);
   else
      error('globalorder:noexactsolver', 'Problem does not have an exact solution available');
   end
elseif strcmp(exactsolver, 'ode15s'),
   disp([datestr(clock, 31), ' |-> Starting reference solution...']);
   disp([' Exact solver: ', exactsolver]);
   ts_ref = [min(tspan), max(tspan)];
   ts_ref = [ts_ref(1), sum(ts_ref)/2, ts_ref(2)];
   
   t = cputime;
   [ts, refsols] = ode15s(problem.LplusN, ts_ref, ...
       problem.y0, ode15soptions, problem);
   t = cputime - t;

   refsol = refsols(end,:);
   disp([datestr(clock, 31), ' |-> Reference solution computed in ', ...
       num2str(t, '%0.4g'), ' seconds']);

else
   disp([datestr(clock, 31), ' |-> Starting reference solution...']);
   % Use an expintscheme for the exact solution, but with a timestep
   % much smaller than the smallest timestep used in the global order
   % test.
   hmin = min(dt) / dtMult;
   disp([' Exact solver: ', exactsolver, ', h = ', num2str(hmin)]);
   
   t = cputime;
   [ts, refsol] = expglm(problem, tspan, hmin, exactsolver);
   t = cputime - t;
   
   refsol = refsol(end,:);
   disp([datestr(clock, 31), ' |-> Reference solution computed in ', ...
       num2str(t, '%0.4g'), ' seconds']);

end


% Some massaging of the schemes-cell-structure:
relstages = ones([1, numel(schemes)]);
for k = 1:numel(schemes),

   % Make space for the error-results:
   schemes{k}.error = zeros([Nexp, 1]);

   % Make space for timing-results:
   schemes{k}.timing = zeros([numel(dt), 1]);

   % We insert the timesteps used in the scheme structure as well.
   schemes{k}.hs = dt;

   % Grab out relstages if present for the schemes.
   % This is if some schemes need smaller timestep in order
   % to be compared fairly.
   if isfield(schemes{k}, 'relstages'),
      relstages(k) = schemes{k}.relstages;
   end
end

progress = cumsum(diff(tspan) ./ dt);
progress = floor(100 * progress ./ progress(end));

siz = size(refsol);

cput = cputime;

% Do for all supplied timesteps:
for e = 1:Nexp,
   expstring = [datestr(clock, 31),                     ...
                ' |-> Experiment ', num2str(e, '%03d'), ...
                ', dt = ', num2str(dt(e), '%0.5e ...')];
   fprintf('%s', expstring);

   % Do for each chosen scheme:
   for k = 1:numel(schemes),
      cpui = cputime;
      [t, numsol] = expglm(problem, tspan, dt(e) / relstages(k), ...
                               schemes{k}.coffcn);
      schemes{k}.timing(e) = cputime - cpui;
      schemes{k}.error(e) = norm(refsol - reshape(numsol, siz), inf)/norm(refsol, ...
                                                        inf);
   end

   fprintf(' (%3d%% done)\n', progress(e));
end


cput = cputime - cput;
disp([datestr(clock, 31), ' |-> Experiments completed in ', ...
       num2str(cput, '%0.4g'), ' seconds']);
