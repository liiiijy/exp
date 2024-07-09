function u = sinegordon_post(y, problem)
% SINEGORDON_POST - Picks out the coordinate value of the phase space
%                   value in the Hamiltonian sine-Gordon problem.
%                   Transforms from Fourier space.
%
% SYNOPSIS:
%   H = sinegordon_post(y, problem);
% 
% DESCRIPTION:
%   Returns the coordinate value (first half) of the phase space
%   variable y, and applies the inverse Fourier transform.
%
% PARAMETERS:
%   y       - Evaluation point in Fourier space, y = (q, p).
%   problem - Problem dependent parameters.  Should be defined by SINEGORDON.
%
% RETURNS:
%   u       - The coordinate value, transformed from Fourier space.
%
% SEE ALSO:
%   SINEGORDON, SINEGORDON_HAMILT

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.4 $  $Date: 2005/11/09 03:18:29 $

u = ifft(y(1:problem.ND));
