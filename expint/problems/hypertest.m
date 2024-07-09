function problem = hypertest(varargin)
% HYPERTEST - Generates a problem structure for the nonlinear Schrodinger
%        equation,
%
%            i y_t = y_xx - 1/(1+y^2) + Phi
%
% SYNOPSIS:
%   problem = hypertest;
%   problem = hypertest(problemdata);
%   problem = hypertest('f1', v1, ..., 'fn', vn);
%
% PARAMETERS:
%   problemdata - OPTIONAL STRUCT specifying problem parameters.  The
%                 following fields are accepted:
%                   ND - Number of Fourier modes.  Must be power of 2.
%                        Default value: ND = 256;
%                   IC - Initial condition.  See source code for
%                        applicable names.
%                        Default value: IC = 'smooth';
%
%   'f1', ..., vn -
%                 OPTIONAL field/value pairs similar to STRUCT
%                 constructor.  Valid field names are the same as for
%                 `problemdata'.
%
% RETURNS:
%   problem  - Structure containing problem specific parameters of the HYPERTEST
%         equation.  `problem' has at least the following fields:
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
% $Revision: 1.2 $  $Date: 2005/10/13 13:37:40 $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Handle default and user choices:
user_val = struct([]); 

default_val = struct('ND',   200, ...
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
problem.L = -1i * L ./ (problem.spatiallength / problem.ND)^2;

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
      problem.exactsolution = 'hypertest_exact';
   otherwise
      if strmatch('reg', lower(IC)),
         [problem.y0, regstr] = regular_function(IC, ND);
         dofft = false;
         problem.ICname = ['Reg', regstr];
         problem.ICnametex = ['\mathrm{Reg', regstr, '}'];
      else
         error('hypertest: Unknown IC supplied');
      end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 4) Setting up the problem parts

% The nonlinear function
problem.N = 'hypertest_N';

% Postprocessing of data is the inverse fourier transform
problem.postprocessing = 'trivial_post';

% The function Lu + N(u,t), needed for ode15s.
problem.LplusN = 'hypertest_LplusN';

problem.problemname = ['Hyperbolic test, ND=', int2str(problem.ND), ...
                       ', IC: ', problem.ICname];
	   
problem.problemnametex = ['Hyperbolic test, $\mathrm{ND}=', int2str(problem.ND), ...
                          '$, IC: $', problem.ICnametex, '$' ];
	   
