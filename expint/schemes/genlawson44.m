function [u, v, a, b, c] = genlawson44(z, problem)
% GENLAWSON44 - Generalized Lawson scheme of stiff order 4.
%
% SYNOPSIS:
%   [u, v, a, b, c] = genlawson44(z, problem);
%
% PARAMETERS:
%   z       - Evaluation point (h * linear operator)
%   problem - Problem dependent parameter structure.
%
% RETURNS:
%   u   - The u_ij(z) coefficient functions as a 4-by-4 cell array.
%   a   - The a_ij(z) coefficient functions as a 4-by-4 cell array.
%   v   - The v_ij(z) coefficient function  as a 4-by-4 cell array.
%   b   - The b_ij(z) coefficient functions as a 4-by-4 cell array.
%   c   - The quadrature nodes as a 4-by-1 DOUBLE array.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.9 $  $Date: 2005/10/10 07:22:12 $

[one, ez2, ez, z] = oneez2(z);

[phi_1, phi_2, phi_3, phi_4]     = phipade(z, 4);
[phi_12, phi_22, phi_32, phi_42] = phipade(1/2*z, 4);

u = { one, [], [], []; ...
      ez2, -3/4*phi_22 - 5/8*phi_32 - 3/16*phi_42, ...
      3/8*phi_22 + 1/2*phi_32 + 3/16*phi_42, ...
     -1/12*phi_22 - 1/8*phi_32 - 1/16*phi_42; ...
      ez2,  -3/4*phi_22 - 5/8*phi_32 - 3/16*phi_42 + 35/32*one, ...
      3/8*phi_22 + 1/2*phi_32 + 3/16*phi_42 - 21/32*one, ...
     -1/12*phi_22 - 1/8*phi_32 - 1/16*phi_42 + 5/32*one; ...
      ez,  -3*phi_2 - 5*phi_3 - 3*phi_4 + 35/16*ez2, ...
      3/2*phi_2 + 4*phi_3 + 3*phi_4 - 21/16*ez2, ...
     -1/3*phi_2 - phi_3 - phi_4 + 5/16*ez2 };

a = { [], [], [], []; ...
      1/2*phi_12 + 11/24*phi_22 + 1/4*phi_32 + 1/16*phi_42, [], [], []; ...
      1/2*phi_12 + 11/24*phi_22 + 1/4*phi_32 + 1/16*phi_42 - 35/32*one, ...
      1/2*one, [], []; ...
      phi_1 + 11/6*phi_2 + 2*phi_3 + phi_4 - 35/16*ez2, [], ez2, [] };

v = { ez,-3*phi_2 - 5*phi_3 - 3*phi_4 + 35/24*ez2 + one, ...
      3/2*phi_2 + 4*phi_3 + 3*phi_4 - 7/8*ez2 - 2/3*one, ...
     -1/3*phi_2 - 1*phi_3 - 1*phi_4 + 5/24*ez2 + 1/6*one; ...
      [], [], [], []; ...
      [], one, [], []; ...
      [], [], one, [] };

b = { phi_1 + 11/6*phi_2 + 2*phi_3 + phi_4 - 35/24*ez2 - 2/3*one, ...
      1/3*ez2, 1/3*ez2, 1/6*one; ...
      one, [], [], []; ...
      [], [], [], []; ...
      [], [], [], [] };

c = [ 0; 1/2; 1/2; 1 ];
