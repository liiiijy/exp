function problem = hochostint(varargin)
% HOCHOSTINT - Generates a problem structure for the Hochbruck-Ostermann
%              integro-differential equation.
%
%       u_t  = u_xx + \int_0^1 u(x,t) dx + \Phi(x,t)
%
% SYNOPSIS:
%   problem = hochostint;
%   problem = hochostint(problemdata);
%   problem = hochostint('f1', v1, ..., 'fn', vn);
%
% PARAMETERS:
%   problemdata - OPTIONAL STRUCT specifying problem parameters.  The
%                 following fields are accepted:
%                   ND - Number of subintervals (physical resolution).
%                        Default value: Same as corresponding parameter
%                        of HOCHOST.
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
%   problem  - Structure containing problem specific parameters of the
%              HOCHOSTINT equation.  `problem' has at least the following
%              fields:
%
%            x  - Vector of all points in physical space.
%            y0 - Initial condition evaluated for every x(i).
%            L  - gODE linear operator.
%            N  - Name of function evaluating ODE non-linear operator.
%            problemname -
%                 String with relevant details of specific problem. 
%                 May be used in plot titles and filenames.
%

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.7 $  $Date: 2005/10/13 13:33:30 $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% This problem is just a modification of the Hochbruck-Ostermann-problem,
% so we reuse that:
problem = hochost(varargin{:});

% We only change the nonlinear function:
problem.N = 'hochostint_N';

% And modifiy the problemname:
warning off
problem.problemname = regexprep(problem.problemname, '(bolic)', ...
                                '$1Int', 'tokenize');
warning on
% Warning is turned off because the above issues a warning in
% Matlab R14 where 'tokenize' is default, while we would like to 
% have this option present here for R13.