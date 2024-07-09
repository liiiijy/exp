function problem = hochost(varargin)
% HOCHOST - Problem specific parameters of the Hochbruck--Ostermann problem.
%           This is the problem in equation 6.1 of Hochbruck and Ostermann's
%           2005 paper: "Explicit exponential Runge-Kutta methods for 
%           semilinear parabolic problems", SIAM J. on Num. Anal.
%
%             y_t = y_xx + 1/(1+u^2) + Phi
%
% SYNOPSIS:
%   problem = hochost;
%   problem = hochost(problemdata);
%   problem = hochost('f1', v1, ..., 'fn', vn);
%
% PARAMETERS:
%   problemdata - OPTIONAL STRUCT specifying problem parameters.  The
%                 following fields are accepted:
%                   ND - Number of subintervals (physical resolution).
%                        Default value: ND = 200;
%                   IC - Initial condition.  See source code for
%                        applicable names.
%                        Default value: IC = 'exactic';
%
%   'f1', ..., vn -
%                 OPTIONAL field/value pairs similar to STRUCT
%                 constructor.  Valid field names are the same as for
%                 `problemdata'.
%
% RETURNS:
%   problem  - Structure containing problem specific parameters of the HOCHOST
%              equation.  `problem' has at least the following fields:
%
%            x  - Vector of all points in physical space.
%            y0 - Initial condition evaluated for every x(i).
%            L  - ODE linear operator.
%            N  - Name of function evaluating ODE non-linear operator.
%            problemname -
%                 String with relevant details of specific problem.
%                 May be used in plot titles and filenames.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.8 $  $Date: 2005/10/13 13:37:40 $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Handle default and user choices:
user_val = struct([]);

default_val = struct('ND', 200, ...
                     'IC', 'exactic');

if nargin > 0,
   user_val = makestruct(varargin{:});
end
s = mergestructs(default_val, user_val);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  2) Discretisation and equation specific parameters
problem.ND = s.ND;

% Only the ND-1 inner points here.
problem.spatiallength = 1; % x \in [0,1] in space.
problem.x = ((1:problem.ND-1) ./ problem.ND)';

e = ones(problem.ND-1, 1);
L = spdiags([e -2*e e], -1:1, problem.ND-1, problem.ND-1);
problem.L = L ./ (problem.spatiallength / problem.ND)^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Initial condition
IC = s.IC;
switch lower(IC);
   case {'hat'}
      problem.y0 = abs(problem.x-1/2) - 1/2;
      problem.ICname = 'hat';
      problem.ICnametex = '\|x-1/2\|-1/2';
   case {'smooth'}
      problem.y0 = exp(-10.*sin(problem.x/2).^2);
      problem.ICname = 'exp(sin(2*pi*x))';
      problem.ICnametex = '\exp(\sin(2\pi x))';
   case {'exactic'}
      problem.y0 = problem.x.*(1-problem.x);
      problem.ICname = 'x(1-x)';
      problem.ICnametex = 'x(1-x)';
      % The exact solution of the PDE  (this also depends on N(t,u))
      problem.exactsolution = 'hochost_exact';
   otherwise
      if strmatch('reg', lower(IC)),
         [problem.y0, regstr] = regular_function(IC, ND);
         dofft = false;
         problem.ICname = ['Reg', regstr];
         problem.ICnametex = ['\mathrm{Reg', regstr, '}'];
      else
         error('hochost:invalidic', 'Unknown IC supplied');
      end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4) Setting up the problem parts

% The nonlinear function
problem.N = 'hochost_N';

% Postprocessing of the data
problem.postprocessing = 'trivial_post';

% The function Lu + N(u,t), needed for ode15s.
problem.LplusN = 'hochost_LplusN';

% Descriptive name to be used in plots or filenames.
problem.problemname = ['Hochbruck-Ostermann, ND=', int2str(problem.ND), ...
                       ', IC: ', problem.ICname];

problem.problemnametex = ['Hochbruck--Ostermann, $\mathit{ND}=', int2str(problem.ND), ...
                          '$, IC: $', problem.ICnametex, '$'];
