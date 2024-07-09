function problem = burgers(varargin)
% BURGERS - Generates a problem structure for the 1D Burgers equation.
%
%             y_t = lambda y_xx - 1/2 (y^2)_x
%
% SYNOPSIS:
%   problem = burgers; 
%   problem = burgers(problemdata); 
%   problem = burgers('f1', v1, ..., 'fn', vn); 
%
% PARAMETERS:
%   problemdata - OPTIONAL STRUCT specifying problem parameters.  The
%                 following fields are accepted:
%                   ND - Number of subinterval (physical resolution).
%                        Default value: ND = 128;
%                   IC - Initial condition.  See source code for
%                        applicable names.
%                        Default value: IC = 'smooth';
%               lambda - Strength of linear operator.
%                        Default value: lambda = 0.03; 
%
%   'f1', ..., vn -
%                 OPTIONAL field/value pairs similar to STRUCT
%                 constructor.  Valid field names are the same as for
%                 `problemdata'.
%
% RETURNS:
%   problem  - Structure containing problem specific parameters of the BURGERS
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
% $Revision: 1.8 $  $Date: 2005/10/12 16:28:01 $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 1) Handle default and user choices:
user_val = struct([]);

default_val = struct('ND', 128, ...
                     'IC', 'smooth', ...
                     'lambda', 0.03);

if nargin > 0, 
  user_val = makestruct(varargin{:}); 
end
s = mergestructs(default_val, user_val);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%  2) Discretisation and equation specific parameters
problem.ND = s.ND;

% Only the N-1 inner points here.
problem.x = 2*pi/problem.ND * (-problem.ND/2 : problem.ND/2 - 1)';

% Equation parameters
problem.lambda = s.lambda;

% Set up L matrix
problem.k = [0:problem.ND/2-1 0 -problem.ND/2+1:-1]';
problem.L = -problem.lambda*problem.k.^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Initial condition
IC = s.IC;
switch lower(IC);
   case {'hat'}
      problem.y0 = pi - abs(problem.x);
      problem.ICname = 'hat';
      problem.ICnametex = '\|x\|';
   case {'smooth'}
      problem.y0 = exp(-10 .* sin(problem.x/2).^2);
      problem.ICname = 'exp(-10sin^2(x/2))';
      problem.ICnametex = '\exp(-10\sin^2(x/2))';
   otherwise
      if strmatch('reg', lower(IC)),
	 [problem.y0, regstr] = regular_function(IC, ND);
	 dofft = false;
	 problem.ICname = ['Reg', regstr];
	 problem.ICnametex = ['\mathrm{Reg', regstr, '}'];
      else
	 error('burgers:invalidic', 'Unknown IC supplied');
      end
end
problem.y0 = fft(problem.y0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 4) Setting up the problem parts

% The nonlinear function
problem.N = 'burgers_N';

% Postprocessing of the data
problem.postprocessing = 'burgers_post';

% The function Lu + N(u,t), needed for ode15s.
problem.LplusN = 'diagLplusN';

% Descriptive name to be used in plots or filenames.
problem.problemname = ['Burgers, ND=', int2str(problem.ND), ...
                       ', IC: ', problem.ICname, ...
                       ', \lambda=', num2str(problem.lambda)];

problem.problemnametex = ['Burgers, $\mathit{ND}=', int2str(problem.ND), ...
                          '$, IC: $', problem.ICnametex, ...
                          '$, $\lambda = ', num2str(problem.lambda), '$' ];
