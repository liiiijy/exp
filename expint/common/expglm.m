function [t, u, varargout] = expglm(problem, tspan, h, cof, varargin) 
% EXPGLM --  Exponential integrator for ODEs of the type
%                u_t = Lu + N(u,t). Integrates up to final time, with
%                chosen scheme. 
%
% SYNOPSIS:
%   [t, u]     = expglm(problem, tspan, h, cof);
%   [t, u]     = expglm(problem, tspan, h, cof, tp);
%   [t, u, uf] = expglm(...);
%
% DESCRIPTION:
%   Implements a generic exponential integrator for the ODE system
%
%           u_t = Lu + N(u)           (*)
%        u(t_0) = u0
%
%
% PARAMETERS:
%   problem -
%            Structure describing the problem. Contains at least
%            the linear part, the non-linear part and the initial
%            condition at tspan(1).
%   tspan  - Time span across which integration is desired.
%   h      - Constant step size. diff(tspan)/h must be an integer for
%            multistep schemes to work. You will get erraneous
%            results if this is not satisfied.
%   cof    - MATLAB function implementing the coefficient functions
%
%               uu_{ij}(z), vv_{ij}(z), aa_{ij}(z), and, bb_{ij}(z)
%
%            Must be callable as
%
%               [uu, vv, aa, bb] = cof(z);
%
%            whence
%
%               uu  -- s-by-r CELL array such that  u{i,j} == uu_{ij}(z)
%               vv  -- r-by-r CELL array such that  v{i,j} == vv_{ij}(z)
%               aa  -- s-by-s CELL array such that  a{i,j} == aa_{ij}(z)
%               bb  -- r-by-s CELL array such that  b{i,j} == bb_{ij}(z)
%   tp     - Time points at which output is desired.
%            The output will be provided at integer multiples of h,
%            rather than the exact values provided here. The code
%            guarantees that max(tp - t) <= h/2 where t is the
%            returned vector of timepoints.
%            OPTIONAL.  
%            Default: Provide output at last step only.
%
% RETURNS:
%   t      - Vector of time points at which numerical solution is
%            computed.
%   u      - Vector of numerical solutions at points
%            specified by tp, or at final step. u(n,:) \approx u(t(n)).
%   uf     - Processed solutions (eg. for plotting).
%            OPTIONAL.
%            The problem structure needs a field set if this is wanted.
%            Cell-array with one entry pr. wanted timepoint.
% LIMITATIONS:
%
%   - Only constant step sizes are supported.
%

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.6 $ $Date: 2005/10/27 19:40:57 $


u0 = problem.y0;

% Call the scheme's coefficient functions
[uu, vv, aa, bb, c] = feval(cof, h .* problem.L, problem);

% r is the number of steps in the scheme. r=1 for
% exponential-Runge-Kutta-schemes.
r = size(bb, 1);

% s is the number of internal stages in the last step.
s = size(aa, 1);

% Determine time step sequence.  Explicitly support a short last step,
% but this will only work for one-step schemes. Excpect erraneous
% results if you supply non-integer multiples of h and use multistep
% methods.
tspan = [min(tspan), max(tspan)];
NSteps = round(diff(tspan) ./ h);
lastshort = false; % Boolean flag saying whether the last step is shorter than
                   % the other.

if ~floatequals(tspan(1) + NSteps*h, tspan(2), 1.0e-10, 5.0e-11),
   lastshort = true;

   last = mod(diff(tspan), h);
   if tspan(1) + h * NSteps > tspan(2),
      NSteps = NSteps - 1;
   end
   % Multisteps do not support taking a short last step, so
   % make sure for multistep methods that that all steps are constant
   if r > 1,
      error('expglm:multistepconstanth', ...
	  ['For multistep methods all stepsizes must be constant.', ...
             ' Change supplied tspan or h']);
   end
end

if nargin > 4,
   % tp is the vector containing timepoints at which 
   % the solution is to be returned.
   tp = varargin{1};
   
   % Check that inputted timepoints are inside tspan:
   if ~(min(tp) >= (tspan(1)) && max(tp) <= tspan(end)),
      error(['expglm:invalidtimepoints:Supplied timepoints are not' ...
	     ' inside tspan']);
   end
   
   % The `round()' here leads to a possible deviation less than 
   % h/2 in the returned timepoints compared to what you want to
   % be returned.
   tpint = round((tp - tspan(1)) ./ h);

   % Note: The actual returned vector of timepoints is built
   % dynamically in the time-loop beneath, based on tpint.

   if ~floatequals(tpint .* h + tspan(1), tp),
      warning('expglm:incompatibletimesteps', ['Supplied' ...
	     ' timestep does not divide all wanted timepoints']);
   end

   % If wanted timepoints are supplied, then we do not need
   % lastshort-functionality, as we are rounding returned
   % timepoints anyway.
   lastshort = false;
else
   tp = tspan(end);
   tpint = NSteps;
end

% Preallocate necessary resources...
u = complex(zeros([length(u0), length(tp)]));
t = zeros([1, length(tp)]);
stage = complex(zeros([length(u0), s]));
out = complex(zeros([length(u0), r-1]));

n = 0;
tn = tspan(1);

% Starting procedure, ONLY done for true multistep methods.
if r > 1,
  
  if isfield(problem, 'exactsolution')
    for l = 1:r-1,
      out(:,r-l) = h .* feval(problem.N, u0, tn, problem);

      tn = tn + h;
      n = n + 1;

      u0 = feval(problem.exactsolution, problem.x, tn, problem);

      t(1) = tn;
      u(:,1) = u0;
    end
  else
    [uuu, vvv, aaa, bbb, cc] = feval('hochost4', h .* problem.L, problem);
    ss = size(aaa, 1);
    
    for l = 1:r-1,
      out(:,r-l) = h .* feval(problem.N, u0, tn, problem);
      
      for in = 1:ss,
        U = uuu{in} * u0;
        
        for j = 1:in-1,
          if ~isempty(aaa{in,j}), U = U + aaa{in,j} * stage(:,j); end
        end
        
        stage(:,in) = h .* feval(problem.N, U, tn + cc(in)*h,  problem);
      end
      tn = tn + h;
      n = n + 1;

      u0 = vvv{1} * u0;

      for ex = 1:ss,
        if ~isempty(bbb{ex}), u0 = u0 + bbb{ex} * stage(:,ex); end
      end

      t(1) = tn;
      u(:,1) = u0;
    end
  end
end

% This is the main stepping routine once the starting procedure
% has been constructed.

out = [u0, out];

for k = 1:length(tp),
   while n < tpint(k),

      % Internal stages
      for in = 1:s,
         U = uu{in,1} * out(:,1);

         for j = 2:r
            if ~isempty(uu{in,j}), U = U + uu{in,j} * out(:,j); end
         end

         for j = 1:in-1,
            if ~isempty(aa{in,j}), U = U + aa{in,j} * stage(:,j); end
         end

         stage(:,in) = h .* feval(problem.N, U, tn + c(in)*h, problem);
      end

      tn = tn + h;
      n = n + 1;

      % Output (GLM-term) approximations:
      for ex = 1:r,
         V = complex(zeros(size(u0)));

         for j = 1:r
            if ~isempty(vv{ex,j}),  V = V + vv{ex,j} * out(:,j); end
         end

         for j = 1:s,
            if ~isempty(bb{ex,j}), V = V + bb{ex,j} * stage(:,j); end
         end

         outer(:,ex) = V;
      end

      out = outer;

      % If the problem-field
      if isfield(problem, 'outputfcn') && ischar(problem.outputfcn),
         feval(problem.outputfcn, out(:,1), tn, problem);
      end
   end

   t(k) = tn;
   u(:,k) = out(:,1);
end

% Note: This code is an (almost) exact duplicate of the above sequence
% but is nontheless repeated here to benefit the most from JIT'ing (R13
% and later).
%
% In a less critical path this would ideally be refactored into a common
% single step function.

if lastshort,
   [uu, vv, aa, bb, cc] = feval(cof, last .* problem.L, problem);

   for in = 1:s,
      U = uu{in,1} * out(:,1);

      for j = 2:r
         if ~isempty(uu{in,j}), U = U + uu{in,j} * out(:,j); end
      end

      for j = 1:in-1,
         if ~isempty(aa{in,j}), U = U + aa{in,j} * stage(:,j); end
      end

      stage(:,in) = last .* feval(problem.N, U, ...
                                  tn + c(in)*last, problem);
   end

   u(:,end) = vv{ex,1} * out(:,1);

   for j = 2:r,
      if ~isempty(vv{ex,j}),
          u(:,end) = u(:,end) + vv{ex,j} * out(:,j);
      end
   end

   for j = 1:s,
      if ~isempty(bb{1,j}),
         u(:,end) = u(:,end) + bb{1,j} * stage(:,j);
      end
   end

   t(end) = tspan(2);
end

% Return the u-vector in the same format as odesuite does:

u = u.';

% This is for making a "processed" solution, for example
% the absolute value for complex solution, making something
% suitable for plotting.

if isfield(problem, 'postprocessing') && nargout > 2,
   for k = 1:numel(tp),
      uf{k} = feval(problem.postprocessing, u(k,:), problem);
   end

   varargout{1} = uf;
end
