function problem = kursiv(varargin)
% KURSIV - Generates a problem structure for the Kuramoto-Sivashinsky 
%          equation.
%
%             y_t = - y_xx - y_xxxx - yy_x
%
% SYNOPSIS:
%   problem = kursiv; 
%   problem = kursiv(problemdata); 
%   problem = kursiv('f1', v1, ..., 'fn', vn); 
%
% PARAMETERS:
%   problemdata - OPTIONAL STRUCT specifying problem parameters.  The
%                 following fields are accepted:
%
%                   ND - Number of subinterval (physical resolution).
%                        Default value: ND = 128;
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
%   problem  - Structure containing problem specific parameters of the KURSIV
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
user_val = struct([]);  % Use default values by default ;-)

default_val = struct(...
    'ND', 128, ...
    'IC', 'smooth');

if nargin > 0, 
   user_val = makestruct(varargin{:}); 
end
s = mergestructs(default_val, user_val);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Discretisation and equation specific parameters
problem.ND = s.ND;
if bitand(problem.ND, problem.ND-1), 
  error('kursiv:ND', 'ND must be a power of two'); end

% Only the N-1 inner points here.
problem.x = 32*pi*(1:problem.ND)'/problem.ND;

% Set up L matrix
problem.k = [0:problem.ND/2-1 0 -problem.ND/2+1:-1]'/16;
problem.L = problem.k.^2 - problem.k.^4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Initial condition
IC = s.IC;
switch lower(IC);
   case {'smooth'}
      problem.y0 = cos(problem.x/16).*(1+sin(problem.x/16));
      problem.ICname = 'cos(x/16)(1+sin(x/16))';
      problem.ICnametex = '\cos(x/16)(1+\sin(x/16))';
   case {'gauss'}
      problem.y0 = exp(-problem.x.^2 / 2);
      problem.ICname = 'gauss'; 
      problem.ICnametex = '\mathrm{exp}(-x^2/2)';
   case {'hat'}
      problem.y0 = 32*pi - abs(problem.x);
      problem.ICname = 'hat';
      problem.ICnametex = '\text{hat}';
   otherwise
      if strmatch('reg', lower(IC))
	 [problem.y0, regstr] = randdecayfcn(IC, problem.ND);
	 problem.ICname = ['Reg', regstr];
	 problem.ICnametex = ['\mathrm{Reg', regstr, '}'];
      else
	 error('kursiv:invalidic', 'Unknown IC supplied');
      end
end
problem.y0 = fft(problem.y0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 4) Setting up the problem parts

% The nonlinear function
problem.N = 'kursiv_N';

% Postprocessing of the data
problem.postprocessing = 'kursiv_post';

% The function Lu + N(u,t), needed for ode15s.
problem.LplusN = 'diagLplusN';

problem.problemname = ['Kuramoto-Sivashinsky, ND=', int2str(problem.ND), ...
                       ', IC: ', problem.ICname];

problem.problemnametex = ['Kuramoto--Sivashinsky, ', ...
                          '$\mathrm{ND}=', int2str(problem.ND), ...
                          '$, IC: $', problem.ICnametex, '$'];
