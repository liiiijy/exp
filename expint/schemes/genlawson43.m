function [u, v, a, b, c] = genlawson43(z, problem)
% GENLAWSON43 - Generalized Lawson scheme of stiff order 4.
%
% SYNOPSIS:
%   [u, v, a, b, c] = genlawson43(z, problem);
%
% PARAMETERS:
%   z       - Evaluation point (h * linear operator)
%   problem - Problem dependent parameter structure.
%
% RETURNS:
%   u   - The u_ij(z) coefficient functions as a 4-by-3 cell array.
%   a   - The a_ij(z) coefficient functions as a 4-by-4 cell array.
%   v   - The v_ij(z) coefficient function  as a 3-by-3 cell array.
%   b   - The b_ij(z) coefficient functions as a 3-by-4 cell array.
%   c   - The quadrature nodes as a 4-by-1 DOUBLE array.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.8 $ $Date: 2005/05/11 11:00:14 $

[one, ez2, ez, z] = oneez2(z);

[phi_1, phi_2, phi_3]    = phipade(z, 3);
[phi_12, phi_22, phi_32] = phipade(1/2*z, 3);

u = { one,                             [],                             []; ...
      ez2,         -1/2*phi_22-1/4*phi_32,          1/8*phi_22+1/8*phi_32; ...
      ez2, -1/2*phi_22-1/4*phi_32+5/8*one, 1/8*phi_22+1/8*phi_32-3/16*one; ...
       ez,   -2*phi_2 - 2*phi_3 + 5/4*ez2,    1/2*phi_2 + phi_3 - 3/8*ez2 };

a = {                                              [],      [],   [], []; ...
                 1/2*phi_12 + 3/8*phi_22 + 1/8*phi_32,      [],   [], []; ...
     1/2*phi_12 + 3/8*phi_22 + 1/8*phi_32 - 15/16*one, 1/2*one,   [], []; ...
                 phi_1 + 3/2*phi_2 + phi_3 - 15/8*ez2,      [],  ez2, [] };

v = { ez, -2*phi_2 - 2*phi_3 + 5/6*ez2 + 1/2*one, 1/2*phi_2 + ...
       phi_3 - 1/4*ez2 - 1/6*one; ...
      [],  [], []; ...
      [], one, [] };

b = {phi_1 + 3/2*phi_2 + phi_3 - 5/4*ez2 - 1/2*one, 1/3*ez2, 1/3*ez2, 1/6*one; ...
                                               one,      [],      [],      []; ...
                                                [],      [],      [], [] };

c = [ 0; 1/2; 1/2; 1 ];
