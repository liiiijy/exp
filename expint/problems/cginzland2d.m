function problem = cginzland2d(varargin)
% CGINZLAND2D - Generates a problem structure for the 2D complex
%               Ginzburg-Landau equation.
%
%             y_t = y + (1+ib) nabla^2 y - (1+ia)y |y|^2
%
% SYNOPSIS:
%   problem = cginzland2d;
%   problem = cginzland2d(problemdata);
%   problem = cginzland2d('f1', v1, ..., 'fn', vn);
%
% PARAMETERS:
%   problemdata - OPTIONAL STRUCT specifying problem parameters.  The
%                 following fields are accepted:
%                   ND - Number of Fourier modes in each spatial
%                        direction.  Must be power of 2.
%                        Default value: ND = 128;
%                   IC - Initial condition.  See source code for
%                        applicable names.
%                        Default value: IC = 'gaussianpulses';
%                   a -  Nonlinear constant.
%                        Default value: a = 1.3
%                   b -  constant in linear part
%                        Default value: b = 0;
%
%   'f1', ..., vn -
%                 OPTIONAL field/value pairs similar to STRUCT
%                 constructor.  Valid field names are the same as for
%                 `problemdata'.
%
% RETURNS:
%   problem  - Structure containing problem specific parameters of the
%              CGINZLAND2D equation.  `problem' has at least the
%              following fields:
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
% $Revision: 1.9 $  $Date: 2005/11/01 14:41:24 $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Handle default and user choices:
user_val = struct([]);

default_val = struct('ND', 128,              ...
                     'IC', 'gaussianpulses', ...
                     'alpha', 0,             ...
                     'beta', 1.3,            ...
                     'length', 200);   % spatial domain size (1D)

if nargin > 0,
   user_val = makestruct(varargin{:});
end
s = mergestructs(default_val, user_val);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  2) Discretisation and equation specific parameters
problem.ND     = s.ND;
problem.alpha  = s.alpha;
problem.beta   = s.beta;
problem.length = s.length;

% The diagonal Laplacian in Fourier space:
m = problem.ND / 2;
kxs = fftshift(linspace(-m, m-1, 2*m))' * (2*pi/problem.length);
kxs(m + 1) = 0;         % set Nyquist frequencies to zero
kys = kxs;

% Spatial discretisation points, ND in x- and y-direction.
problem.x = linspace(0, problem.length, problem.ND)';
problem.y = problem.x;

% Matrix array of Fourier modes:
% L is diagonal, make it a column vector of size ND^2
[KX, KY]  = meshgrid(kxs, kys);
nabla2d   = KX.*KX + KY.*KY;            % The 2D nabla^2
problem.L = ones(problem.ND^2, 1) - (1 + 1i*problem.alpha) .* nabla2d(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Initial condition
% Flag determining whether or not we need to perform DFT on the IC
% specified here.
dofft = true;
switch lower(s.IC);
   case { 'gaussianpulses' }
      Lambda = problem.length; % parameter for gauss-curves.
      L2 = Lambda / 2;
      L3 = Lambda / 3;
      f  = -20 / Lambda;

      [X, Y] = meshgrid(problem.x, problem.y);
      u0 = exp(f .* ((X - L3).^2 + (Y - L3).^2)) - ...
           exp(f .* ((X - L2).^2 + (Y - L2).^2)) + ...
           exp(f .* ((X - L2).^2 + (Y - L3).^2));

      % Fourier values in grid:
      ufourier = fft2(u0);

      % column vector representation of IC
      problem.y0 = ufourier(:);
      problem.ICname = 'GaussianPulses';
      problem.ICnametex = '\text{Gaussian pulses})';
   case { 'random' }
      % Fourier space values on spatial grid:
      ufourier = fft2(1.0e-4 * randn([numel(problem.x), numel(problem.y)]));

      problem.y0 = ufourier(:);
      problem.ICname = 'Random';
      problem.ICnametex = '\text{Random}';
   otherwise
      error('cginzland2d:invalidic', 'Unknown IC supplied');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4) Setting up the problem parts

% The nonlinear function
problem.N = 'cginzland2d_N';

% Postprocessing of data is the inverse fourier transform
problem.postprocessing = 'cginzland2d_post';

% The function Lu + N(u,t), needed for ode15s.
problem.LplusN = 'diagLplusN';

% Descriptive name to be used in plots or filenames.
problem.problemname = ['Complex Ginzburg-Landau 2D, ND=',   ...
                                int2str(problem.ND),        ...
                       ', IC: ', problem.ICname,            ...
                       ', \alpha=', num2str(problem.alpha), ...
                       ', \beta=', num2str(problem.beta)];

problem.problemnametex = ['Complex Ginzburg--Landau, $\mathit{ND}=', ...
                                   int2str(problem.ND), '$, IC: $',  ...
                           problem.ICnametex, '$, $\alpha = ',       ...
                           num2str(problem.alpha), '$, $\beta = ',   ...
                           num2str(problem.beta), '$'];
