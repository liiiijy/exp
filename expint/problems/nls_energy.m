function energy = nls_energy(U, pr)
% NLS_ENERGY --- Compute the energy of a solution vector
%
% SYNOPSIS:
%   energy = nls_energy(U, pr);
% 
% DESCRIPTION:
%   Computes the integral
%    
%     H = \int_-pi^pi  | diff(u(x),x) |^2 - lambda/2 * |u(x)|^4 dx
%
%   with a rectangle-rule.
%
% PARAMETERS:
%    U - a solution value in Fourier space.
%    pr - the problem structure. 
%
% RETURNS:
%    energy - approximated value of the integral.
%
% LIMITATIONS:
%
%   - The potential must be zero.

%
% $Id: nls_energy.m,v 1.2 2005/10/12 16:28:01 berland Exp $
%

% Check that the potential is zero:
if sum(pr.v) ~= 0,
   error('nls_energy:invalidpotential', 'Potential must be zero');
end

% Assume input is Fourier-transformed:
psi = abs(ifft(U)); 

% psi is periodic, append for the use of the diff-command:
psiappended= [psi;  psi(1)];

% Differentiate:
h = 2*pi / pr.ND;

psid = h*diff(psiappended);

energy = sum( abs(psid).^2 - pr.lambda/2* abs(psi).^4 );

