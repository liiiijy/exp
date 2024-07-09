function nls_plotfourier(u, t,  problem)
% NLS_PLOTFOURIER -- Output plotting function for the NLS problem.
%
% SYNOPSIS:
%   nls_plotfourier(u, t, problem)
%
% DESCRIPTION:
%   Plots the absolute value of the Fourier coefficients of the
%   (numerical) solution to the NLS problem, together with 
%   two random function of regularity 1 and 10 for comparison.
%
%   This function is called from EXPGLM as a callback functionality.
%
% PARAMETERS:
%   u       - (Numerical) solution to the NLS problem.
%   t       - point in time. Not used.
%   problem - Problem description structure as defined by NLS.
%
% RETURNS:
%   Nothing.
%
% SEE ALSO:
%   NLS, EXPGLM.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.3 $  $Date: 2005/10/22 02:50:13 $

idxs = find(pr.ks>0); % Ignore negative fourier coefficients.

N = feval(pr.N, u, 0, pr); % time should not be zero...
hold off

reg1ex = abs(randdecayfcn(1, pr.ND));
loglog(pr.ks(idxs), reg1ex(idxs), 'g.');
hold on 

% Plot the solution U in blue:
loglog(pr.ks(idxs), abs(u(idxs)), 'b.');


% Plot the N-function in red:
loglog(pr.ks(idxs), abs(N(idxs)), 'r.');

reg10ex = abs(randdecayfcn(10, pr.ND));
loglog(pr.ks(idxs), reg10ex(idxs), 'g.');

legend('reg1', 'solution', 'N', 'reg10');

title(['Time: ', num2str(t)]);

drawnow
pause
