function problem = allencahn(varargin)
% ALLENCAHN - Generates a problem structure for the Allen-Cahn equation.
%
%             y_t = lambda y_xx + y - y^3
%
% SYNOPSIS:
%   problem = allencahn;
%   problem = allencahn(problemdata);
%   problem = allencahn('f1', v1, ..., 'fn', vn); 
%
% PARAMETERS:
%   problemdata - OPTIONAL STRUCT specifying problem parameters.  The
%                 following fields are accepted:
%                   ND - Number of subintervals (physical resolution).
%                        Default value: ND = 50; 
%                   IC - Initial condition.  See source code for
%                        applicable names.
%                        Default value: IC = 'gaussianpulses';
%               lambda - Strength of linear operator.
%                        Default value: lambda = 0.001; 
%
%   'f1', ..., vn -
%                 OPTIONAL field/value pairs similar to STRUCT
%                 constructor.  Valid field names are the same as for
%                 `problemdata'.
%
% RETURNS:
%   problem  - Structure containing problem specific parameters of the
%              ALLEN-CAHN equation.  `problem' has at least the following
%              fields. 
%
%           x  - Vector of all points in physical space.
%           y0 - Initial condition evaluated for every x(i).
%           L  - ODE linear operator.
%           N  - Name of function evaluating ODE non-linear operator.
%           problemname -
%                String with relevant details of specific problem. 
%                May be used in plot titles and filenames.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.8 $  $Date: 2005/10/13 13:37:40 $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 1) Handle default and user choices:
user_val = struct([]);

default_val = struct('ND', 50, ...
                     'lambda', 0.001, ...
                     'IC', 'smooth');

if nargin > 0, 
  user_val = makestruct(varargin{:}); 
end
s = mergestructs(default_val, user_val);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 2) Discretisation and equation specific parameters
problem.ND = s.ND;

% This is to be interpreted as ND-1 inner points!
problem.spatiallength = 2; % x \in [-1,1] in space.

% Only the N-1 inner points here.
problem.x = cos(pi*(0:problem.ND)/problem.ND)';
problem.x = problem.x(2:problem.ND);

% Equation parameters
problem.lambda = s.lambda;

% Set up L matrix
c = [2; ones(problem.ND-1,1); 2].*(-1).^(0:problem.ND)';
X = repmat([1; problem.x; -1], 1, problem.ND+1);
dX = X - X';                  
D  = (c*(1./c)')./(dX+(eye(problem.ND+1)));      % off-diagonal entries
D  = D - diag(sum(D'));                          % diagonal entries
L  = D * D;
problem.L = problem.lambda*L(2:end-1, 2:end-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Initial condition
IC = s.IC;
switch lower(IC);
   case {'smooth'} % use w = u-x to make BCs homogenous
      problem.y0 = 0.53*problem.x + 0.47*sin(-1.5*pi*problem.x) - problem.x;  
      problem.ICname = '0.53x+0.47sin(-1.5*pi*x)';
      problem.ICnametex = '0.53x+0.47\sin(-1.5\pi x)';
   otherwise
      error('allencahn:invalidic', 'Unknown IC supplied');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 4) Setting up the problem parts

% The nonlinear function
problem.N = 'allencahn_N';

% Postprocessing of the data
problem.postprocessing = 'allencahn_post';

% The function Lu + N(u,t), needed for ode15s.
problem.LplusN = 'matrixLplusN';

% Descriptive name to be used in plots or filenames.
problem.problemname = ['Allen-Cahn, ND=', int2str(problem.ND), ...
                       ', IC: ', problem.ICname, ...
                       ', \lambda=', num2str(problem.lambda)];

problem.problemnametex = ['Allen--Cahn, $\mathrm{ND}=', int2str(problem.ND), ...
                          '$, IC: $', problem.ICnametex, ...
                          '$, $\lambda = ', num2str(problem.lambda), '$'];
