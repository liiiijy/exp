function problem = kdv(varargin)
% KDV - Generates a problem structure for the KDV problem.
%
%             y_t = - y_xxx - yy_x
%
% SYNOPSIS:
%   problem = kdv; 
%   problem = kdv(problemdata); 
%   problem = kdv('f1', v1, ..., 'fn', vn); 
%
% PARAMETERS:
%   problemdata - OPTIONAL STRUCT specifying problem parameters.  The
%                 following fields are accepted:
%                   ND - Number of Fourier modes in each spatial
%                        direction.  Must be power of 2.
%                        Default value: ND = 128;
%                   IC - Initial condition.  See source code for
%                        applicable names.
%                        Default value: IC = 'smooth';
%               lambda - Initial condition amplitude and domain
%                        stretching factor.
%                        Default value: lambda = 625;
%
%   'f1', ..., vn -
%                 OPTIONAL field/value pairs similar to STRUCT
%                 constructor.  Valid field names are the same as for
%                 `problemdata'.
%
% RETURNS:
%   problem  - Structure containing problem specific parameters of the KDV
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
% $Revision: 1.10 $  $Date: 2005/10/12 16:28:01 $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Handle default and user choices:
user_val = struct([]);

default_val = struct('ND', 128, ...
                     'IC', 'smooth', ...
                     'lambda', 625);

if nargin > 0, 
   user_val = makestruct(varargin{:}); 
end
s = mergestructs(default_val, user_val);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Discretisation and equation specific parameters
problem.ND = s.ND;
problem.lambda = s.lambda;

if bitand(problem.ND, problem.ND-1), 
  error('kdv:ND', 'ND must be a power of two'); end

% Only the N-1 inner points here.
problem.x = 2*pi/problem.ND*(-problem.ND/2 : problem.ND/2 - 1)';

% Set up L matrix
problem.k = [0 : problem.ND/2 - 1, 0, -problem.ND/2+1 : -1]';
problem.L = 1i * problem.k .^ 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Initial condition
IC = s.IC;
switch lower(IC);
   case {'smooth'}
      problem.origy0 = 3*problem.lambda*sech(sqrt(problem.lambda)*problem.x/2).^2;
      problem.y0 = fft(problem.origy0 / norm(problem.origy0));
      problem.ICname = 'sech(x)^2';
      problem.ICnametex = '\mathrm{sech}(x)^2';
   otherwise
      error('kdv:invalidic', 'Unknown IC supplied');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 4) Setting up the problem parts

% The nonlinear function
problem.N = 'kdv_N';

% Postprocessing of the data
problem.postprocessing = 'kdv_post';

% The function Lu + N(u,t), needed for ode15s.
problem.LplusN = 'diagLplusN';

% Descriptive name to be used in plots or filenames.
problem.problemname = ['Korteweg-de Vries, ND=', int2str(problem.ND), ...
                       ', IC: ', problem.ICname, ...
                       ', \lambda=', num2str(problem.lambda)];

problem.problemnametex = ['Korteweg--de Vries, $\mathit{ND}=', ...
                          int2str(problem.ND), ...
                          '$, IC: $', problem.ICnametex, ...
                          '$, $\lambda=', num2str(problem.lambda), '$'];
