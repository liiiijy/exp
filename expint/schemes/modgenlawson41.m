function [u, v, a, b, c] = modgenlawson41(z, problem)
% MODGENLAWSON41 - Modified generalized Lawson scheme of order 2.
%
% SYNOPSIS:
%   [u, v, a, b, c] = modgenlawson41(z, problem);
%
% PARAMETERS:
%   z       - Evaluation point (h * linear operator)
%   problem - Problem dependent parameter structure.
%
% RETURNS:
%   u   - The u_ij(z) coefficient functions as a 4-by-1 cell array.
%   a   - The a_ij(z) coefficient functions as a 4-by-4 cell array.
%   v   - The v_ij(z) coefficient function  as a 1-by-1 cell array.
%   b   - The b_ij(z) coefficient functions as a 1-by-4 cell array.
%   c   - The quadrature nodes as a 4-by-1 DOUBLE array.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.7 $  $Date: 2005/05/11 08:30:29 $

[one, ez2, ez, z] = oneez2(z);

[phi_1, phi_2] = phipade(z, 2);
[phi_12]       = phipade(1/2*z, 1);

u = {one; ez2; ez2; ez};

a = {                  [],      [],   [],  []; ...
               1/2*phi_12,      [],   [],  []; ...
     1/2*phi_12 - 1/2*one, 1/2*one,   [],  []; ...
              phi_1 - ez2,      [],  ez2,  [] };

v = { ez };

b = { phi_1 - phi_2 - 1/3*ez2, 1/3*ez2, 1/3*ez2, phi_2 - 1/3*ez2 };

c = [ 0; 1/2; 1/2; 1 ];
