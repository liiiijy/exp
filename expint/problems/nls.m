function problem = nls(varargin)
% NLS - Generates a problem structure for the nonlinear Schrodinger
%       equation,
%
%            i psi_t = - psi_xx + (V(x) + lambda |psi|^2) psi
%
% SYNOPSIS:
%   problem = nls;
%   problem = nls(problemdata);
%   problem = nls('f1', v1, ..., 'fn', vn);
%
% PARAMETERS:
%   problemdata - OPTIONAL STRUCT specifying problem parameters.  The
%                 following fields are accepted:
%                   ND - Number of Fourier modes.  Must be power of 2.
%                        Default value: ND = 256;
%                   IC - Initial condition.  See source code for
%                        applicable names.
%                        Default value: IC = 'smooth';
%                   Potential -
%                        Potential function.  See source code for
%                        applicable names.
%                        Default value: Potential = 'smooth';
%                   lambda -
%                        Nonlinear constant.
%                        Default value: lambda = 0; (linear)
%
%   'f1', ..., vn -
%                 OPTIONAL field/value pairs similar to STRUCT
%                 constructor.  Valid field names are the same as for
%                 `problemdata'.
%
% RETURNS:
%   problem  - Structure containing problem specific parameters of the NLS
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
% $Revision: 1.19 $  $Date: 2005/11/09 03:11:49 $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File organisation:
%   1) Handle default and user choices
%   2) Discretisation and equation specific parameters
%   3) Initial condition
%   4) Potential function
%   5) Structure field specification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Handle default and user choices:
user_val = struct([]); 

default_val = struct( ...
    'ND',   256, ...
    'lambda', 0,  ...
    'IC', 'smooth', ...
    'Potential', 'smooth');

if nargin > 0, 
   user_val = makestruct(varargin{:}); 
end
s = mergestructs(default_val, user_val);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  2) Discretisation and equation specific parameters
ND = s.ND;

if bitand(ND, ND-1), 
   error('nls:ND', 'ND must be a power of two'); 
end

% The diagonal Laplacian in Fourier space:
m  = ND / 2;
ks = fftshift(linspace(-m, m-1, 2*m))';
ks(m + 1) = 0;
L = -1i * ks.^2;

% Equation parameters
lambda = s.lambda;
x   = pi / m * (-m : m-1)'; % The domain has length 2*pi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Initial condition
% Flag determining whether or not we need to perform DFT on the IC
% specified here.  Set to FALSE if RANDDECAYFCN is used.
dofft = true;
IC = s.IC;
switch lower(IC);
   case {'hat'}
      y0 = pi - abs(x);
      ICname = 'hat';
      ICnametex = '\mathrm{hat}';
   case {'exp(sin(2x))', 'smooth'}
      y0 = exp(sin(2 * x));
      ICname = 'exp(sin(2x))';
      ICnametex = '\exp(\sin(2x))';
   case {'sinsqr', '1oversinsqr'}
      y0 = 1 ./ (1 + sin(x).^2);
      ICname = '1overSinSqr';
      ICnametex = '1/(1+\sin^2(x))';
   case {'constant'}
      y0 = ones(length(x), 1);
      ICname = 'constant';
      ICnametex = '1';
   case {'constantfourier'} % Constant 1 in Fourier domain.
      y0 = ones(length(x), 1);
      ICname = 'constantfourier';
      ICnametex = '\mathrm{constantfourier}';
      dofft = false;
   case { 'fracreg' }
      y0 = randdecayfcn(1.5, ND);
      dofft = false;
      ICname = 'Reg1p5';
      ICnametex = '\mathrm{Reg}1.5';
   case { 'perturbedplanewave' }
      y0 = 1/2 * (1 + 0.01 * cos( x ));
      ICname = 'PerturbedPlaneWave';
      ICnametex = '\tfrac12(1+\varepsilon \cos(\mu x))';
   otherwise
      if strmatch('reg', lower(IC)),
	 [y0, regstr] = randdecayfcn(IC, ND);
	 dofft = false;
	 ICname = ['Reg', regstr];
	 ICnametex = ['\mathrm{Reg', regstr, '}'];
      else
	 error('nls:invalidid', 'Unknown IC supplied');
      end
end

if dofft,
   y0 = y0 ./ norm(y0);
   y0 = fft(y0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  4) Potential function:
Potential = s.Potential;
switch lower(Potential);
   case {'constant'}
      v = ones([ND, 1]);
      potname = '1';
      potnametex = potname;
   case {'zero'}
      v = zeros([ND, 1]);
      potname = '0';
      potnametex = potname;
   case {'smooth'}
      v = 1./(1 + sin(x) .^ 2);
      potname = '1overSinSqr';
      potnametex = '1/(1+\sin^2(x))';
   case {'hat'}
      v = abs(x);
      potname = 'abs(x)';
      potnametex = '\abs(x)';
   otherwise
      if strmatch('reg', lower(Potential))
	 [v, regstr] = randdecayfcn(Potential, ND);
	 v = ifft(v);
	 potname = ['Reg', regstr];
	 potnametex = ['\mathrm{Reg', regstr, '}'];
      else
	 error('nls:potentialinvalid', 'Unknown potential supplied');
      end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  5) Setting structure fields

problem.lambda = lambda;
problem.ND = ND;
problem.ks = ks; % may be of some use.
problem.x  = x;
problem.ks = ks;
problem.y0 = y0;  problem.ICname        = ICname;
problem.v  = v;   problem.potentialname = potname;
problem.L = L;

% The nonlinear function
problem.N = 'nls_N';

% Postprocessing of data is the inverse fourier transform
problem.postprocessing = 'nls_post';

% The function Lu + N(u,t), needed for ode15s.
problem.LplusN = 'diagLplusN';

% Uncomment to enable incremental output processing
%problem.outputfcn = 'nls_plotenergy';
%problem.outputfcn = 'nls_plotfourier';

problem.problemname = ['Nonlinear Schrödinger, ND=', int2str(ND), ...
                       ', IC: ', ICname, ...
                       ', Pot: ', potname, ...
                       ', \lambda=', num2str(lambda)];
	   
problem.problemnametex = ['Nonlinear Schr\"odinger, $\mathrm{ND}=', int2str(ND), ...
                          '$, IC: $', ICnametex, ...
                          '$, Pot: $', potnametex, ...
                          '$, $\lambda = ', num2str(lambda), '$' ];
	   
