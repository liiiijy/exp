function [u, v, a, b, c] = hochost4(z, problem)
% HOCHOST4 - Coefficient functions for the the explicit expint RK(s=5)
%            scheme on page 19 in hochbruck04eer with stiff order 4.
%
% SYNOPSIS:
%   [u, v, a, b, c] = hochost4(z, problem);
%
% PARAMETERS:
%   z       - Evaluation point (h * linear operator)
%   problem - Problem dependent parameter structure.
%
% RETURNS:
%   u   - The u_ij(z) coefficient functions as a 5-by-1 cell array.
%   a   - The a_ij(z) coefficient functions as a 5-by-5 cell array.
%   v   - The v_ij(z) coefficient function  as a 1-by-1 cell array.
%   b   - The b_ij(z) coefficient functions as a 1-by-5 cell array.
%   c   - The quadrature nodes as a 5-by-1 DOUBLE array.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.10 $  $Date: 2005/10/10 07:22:12 $

[one, ez2, ez, z] = oneez2(z);

[phi_1, phi_2, phi_3]    = phipade(z, 3, 13);
[phi_12, phi_22, phi_32] = phipade(1/2*z, 3, 13);

u = {one; ez2; ez2; ez; ez2};

a_52 = 1/2*phi_22 - phi_3 + 1/4*phi_2 - 1/2*phi_32;
a_54 = 1/4*phi_22 - a_52;

a = {                      [],     [],    [],    [],  []; ...
                   1/2*phi_12,     [],    [],    [],  []; ...
          1/2*phi_12 - phi_22, phi_22,    [],    [],  []; ...
              phi_1 - 2*phi_2,  phi_2, phi_2,    [],  []; ...
     1/2*phi_12-2*a_52 - a_54,   a_52,  a_52,  a_54,  [] };

v = { ez };

b = { phi_1 - 3*phi_2 + 4*phi_3, [], [], -phi_2 + 4*phi_3, 4*phi_2 - 8*phi_3 };

c = [ 0; 1/2; 1/2; 1; 1/2 ];
