function problem = sinegordon(varargin);
% SINEGORDON - Generates a problem structure for the sine-Gordon
%              equation, spectral discretization in space.
%
%          u_tt = u_xx + sin(u)
%
% The equation is rewritten as a system of first order 
% equations, using z = (u, v),  v = d/dt u. The system is then
% 
%      .       v
%      z  =    u_xx + sin(u)
%
% splitted so that N_1(z) = 0, N_2(z) = sin(u).
%
%
% SYNOPSIS:
%   problem = sinegordon;
%   problem = sinegordon(problemdata);
%   problem = sinegordon('f1', v1, ..., 'fn', vn);
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
%   problem  - Structure containing problem specific parameters of the SINEGORDON
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
% $Revision: 1.5 $  $Date: 2005/11/09 03:18:29 $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Handle default and user choices.
user_val = struct([]); 

default_val = struct( ...
    'ND',  128, ...
    'IC', 'stationarysoliton', ...
    'Length', 40 ...
	);
  

if nargin > 0, 
   user_val = makestruct(varargin{:}); 
end
s = mergestructs(default_val, user_val);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  2) Discretization and equation specific parameters
ND = s.ND;

m  = ND / 2;
dx = s.Length/s.ND;

% Symmetric interval.
x = dx * [-m : (m-1) ]';

problem.ksorig = fftshift(linspace(-m, m-1, 2*m))';
problem.ks = problem.ksorig; % We need to save this ks-vector.
problem.ks(m + 1) = 0; % set nyquist frequency to zero

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Initial condition
switch(lower(s.IC));
   case {'stationarysoliton'} %%% Length should be fairly big, e.g. 40 to
                      % avoid boundary influencing the soliton.
      u0 = zeros(size(x));
      v0 = 4*sech(x);
      y0 = [fft(u0) ; fft(v0)];
   case {'unstable'} % A solution on a homoclinic orbit.
      % Length MUST be 2*sqrt(2)*pi.
      amp = 0.1;
      mu = 2*pi/s.Length;
      u0 = pi + amp*cos(mu*x);
      v0 = zeros(size(u0));
      y0 = [fft(u0) ; fft(v0)];
   otherwise
      error('unknown IC supplied');
end

if (~exist('ICname')), ICname = s.IC; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4) Setting structure fields.

% The linear operator is
%       0   I
%  L =  k^2 0
% and is of size 2ND*2ND.
problem.L                  = spalloc(2*ND, 2*ND, 2*ND);
problem.L(1:ND, ND+1:2*ND) = speye(ND);
problem.L(ND+1:2*ND, 1:ND) = diag((2*pi*i/s.Length * problem.ksorig).^2);
% Note ksorig is used here, not ks.

problem.N      = 'sinegordon_N';
problem.LplusN = 'matrixLplusN';

% The inverse Fourier transform on the first half of the solution vector:
problem.postprocessing = 'sinegordon_post';

problem.dx     = dx;
problem.ND     = ND;
problem.u0     = u0;
problem.y0     = y0;
problem.v0     = v0;
problem.ICname = ICname;
problem.x      = x;      
problem.Length = s.Length;

problem.problemname = ['sine-Gordon, ND=', int2str(ND), ...
                       ', IC: ', ICname, ...
                       ];







