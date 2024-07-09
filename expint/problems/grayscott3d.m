function problem = grayscott3d(varargin)
% GRAYSCOTT3D - Generates a problem structure for the 3D Gray--Scott equation.
%
%             u_t = Du nabla^2 u - uv^2 + (1-u)a
%             v_t = Dv nabla^2 v + uv^2 - (a+b)v
%
% SYNOPSIS:
%   problem = grayscott3d;
%   problem = grayscott3d(problemdata);
%   problem = grayscott3d('f1', v1, ..., 'fn', vn);
%
% PARAMETERS:
%   problemdata - OPTIONAL STRUCT specifying problem parameters.  The
%                 following fields are accepted:
%                   ND - Number of Fourier modes in each spatial
%                        direction.  Must be even.
%                        Default value: ND = 32;
%                   IC - Initial condition.  See source code for
%                        applicable names.
%                        Default value: IC = 'smooth';
%                   Du - Diffusion parameter for `u' species.
%                        Default value: Du = 2.0e-5; 
%                   Dv - Diffusion parameter for `v' species.
%                        Default value: Du = 1.0e-5; 
%                   a -  Nonlinear constant.
%                        Default value: a = 0.035
%                   b -  Nonlinear constant.
%                        Default value: b = 0.065;
%
%   'f1', ..., vn -
%                 OPTIONAL field/value pairs similar to STRUCT
%                 constructor.  Valid field names are the same as for
%                 `problemdata'.
%
% RETURNS:
%   problem  - Structure containing problem specific parameters of the
%              equation.  `problem' has at least the following fields:
%
%            x  - Vector of all points in physical space.
%            y0 - Initial condition as a vector
%            L  - ODE linear operator.
%            N  - Name of function evaluating ODE non-linear operator.
%            problemname -
%                 String with relevant details of specific problem.
%                 May be used in plot titles and filenames.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.8 $  $Date: 2005/11/01 14:41:25 $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Handle default and user choices:
user_val = struct([]);

default_val = struct('ND', 32,       ...
                     'IC', 'smooth', ...
                     'Du', 2.0e-5,   ...
                     'Dv', 1.0e-5,   ...
                     'alpha', 0.035, ...
                     'beta', 0.065,  ...
                     'length', 0.5);    % Spatial domain size (1D).

if nargin > 0,
   user_val = makestruct(varargin{:});
end
s = mergestructs(default_val, user_val);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  2) Discretisation and equation specific parameters
problem.ND     = s.ND;
problem.Du     = s.Du;
problem.Dv     = s.Dv;
problem.alpha  = s.alpha;
problem.beta   = s.beta;
problem.length = s.length;

if mod(problem.ND, 2) ~= 0,
   error('grayscott3d:ND', 'ND must be even');
end

% The diagonal Laplacian in Fourier space:
m   = problem.ND / 2;
kxs = fftshift(linspace(-m, m-1, 2*m))' * (2*pi/problem.length);
kxs(m + 1) = 0;         % Set Nyquist frequencies to zero.
kys = kxs;
kzs = kxs;

% Spatial discretisation points, ND in x- and y-direction.
problem.x = linspace(0, problem.length, problem.ND)';
problem.y = problem.x;
problem.z = problem.x;

% Matrix array of Fourier modes:
% L is diagonal, make it a column vector of size 2*ND^3
[KX, KY, KZ] = meshgrid(kxs, kys, kzs);
nabla3d      = -(KX.*KX + KY.*KY + KZ.*KZ);     % 3D nabla^2
problem.L    = [problem.Du .* nabla3d(:); problem.Dv .* nabla3d(:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Initial condition
% Flag determining whether or not we need to perform DFT on the IC
% specified here.
dofft = true;
switch lower(s.IC);
   case { 'smooth' }
      Lambda = problem.length;  % parameter for gauss-curves.
      L2 = Lambda / 2;
      f  = - 150;

      [X, Y, Z] = meshgrid(problem.x, problem.y, problem.z);
      v0 = exp(f .* ((X - L2).^2 + 10.*(Y - L2).^2 + 5.*(Z - L2).^2));
      u0 = 1 - v0;

      % Fourier values in grid:
      ufourier = fftn(u0);
      vfourier = fftn(v0);

      % column vector representation of concatenated IC values
      problem.y0 = [ufourier(:); vfourier(:)];
      problem.ICname = 'Smooth';
      problem.ICnametex = '\text{Smooth})';
   case { 'random' }
      % Physical values on spatial grid:
      u0 = 1.0e-4*randn([numel(problem.x) numel(problem.y) numel(problem.z)]);
      v0 = 1.0e-4*randn([numel(problem.x) numel(problem.y) numel(problem.z)]);

      ufourier = fftn(u0);
      vfourier = fftn(v0);

      % column vector representation of concatenated IC values
      problem.y0 = [ufourier(:); vfourier(:)];
      problem.ICname = 'Random';
      problem.ICnametex = '\text{Random}';
   otherwise
      error('grayscott3d:invalidic', 'Unknown IC supplied');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4) Setting up the problem parts

% The nonlinear function
problem.N = 'grayscott3d_N';

% Postprocessing of data is the inverse fourier transform
problem.postprocessing = 'grayscott3d_post';

% The function Lu + N(u,t), needed for ode15s.
problem.LplusN = 'diagLplusN';

% Descriptive name to be used in plots or filenames.
problem.problemname = ['Gray-Scott3d, ND=',                 ...
                                int2str(problem.ND),        ...
                       ', IC: ', problem.ICname,            ...
                       ', D_u=', num2str(problem.Du),       ...
                       ', D_v=', num2str(problem.Dv),       ...
                       ', \alpha=', num2str(problem.alpha), ...
                       ', \beta=', num2str(problem.beta)];

problem.problemnametex = ['Gray-Scott3d, $\mathit{ND}=',           ...
                                int2str(problem.ND),               ...
                          '$, IC: $', problem.ICnametex,           ...
                          '$, $D_u = ', num2str(problem.Du),       ...
                          '$, $D_v = ', num2str(problem.Dv),       ...
                          '$, $\alpha = ', num2str(problem.alpha), ...
                          '$, $\beta = ', num2str(problem.beta), '$'];
