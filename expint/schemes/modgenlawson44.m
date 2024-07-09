function [u, v, a, b, c] = modgenlawson44(z, problem)
% MODGENLAWSON44 - Modified generalized Lawson scheme of stiff order 5.
%
% SYNOPSIS:
%   [u, v, a, b, c] = modgenlawson44(z, problem);
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
% $Revision: 1.8 $ $Date: 2005/05/11 11:00:14 $

[one, ez2, ez, z] = oneez2(z);

[phi_1, phi_2, phi_3, phi_4, phi_5] = phipade(z, 5);
[phi_12, phi_22, phi_32, phi_42]    = phipade(1/2*z, 4);

u = { one,  [], [], []; ...
      ez2,  -3/4*phi_22 - 5/8*phi_32 - 3/16*phi_42, ...
      3/8*phi_22 + 1/2*phi_32 + 3/16*phi_42, ...
     -1/12*phi_22 - 1/8*phi_32 - 1/16*phi_42; ...
      ez2,  -3/4*phi_22 - 5/8*phi_32 - 3/16*phi_42 + 35/32*one, ...
      3/8*phi_22 + 1/2*phi_32 + 3/16*phi_42 - 21/32*one, ...
     -1/12*phi_22 - 1/8*phi_32 - 1/16*phi_42 + 5/32*one; ...
      ez,  -3*phi_2 - 5*phi_3 - 3*phi_4 + 35/16*ez2, ...
      3/2*phi_2 + 4*phi_3 + 3*phi_4 - 21/16*ez2, ...
     -1/3*phi_2 - phi_3 - phi_4 + 5/16*ez2 };

a = {                                                   [], [], [], []; ...
      1/2*phi_12 + 11/24*phi_22 + 1/4*phi_32 + 1/16*phi_42, [], [], []; ...
      1/2*phi_12 + 11/24*phi_22 + 1/4*phi_32 + 1/16*phi_42 - 35/32*one, ...
                                                       1/2*one, [], []; ...
      phi_1 + 11/6*phi_2 + 2*phi_3 + phi_4-35/16*ez2,   [], ez2, [] };

v = { ez, -3/2*phi_2 + 1/2*phi_3 + 6*phi_4 + 6*phi_5 + 35/96*ez2, ...
      1/2*phi_2 + 1/3*phi_3 - 3*phi_4 - 4*phi_5 - 7/48*ez2, ...
     -1/12*phi_2 - 1/12*phi_3 + 1/2*phi_4 + phi_5 + 5/192*ez2; ...
      [] , [],  [], []; ...
      [], one,  [], []; ...
      [],  [], one, [] };

b = { phi_1 + 5/6*phi_2 - 5/3*phi_3 - 5*phi_4 - 4*phi_5 - 35/48*ez2, ...
      1/3*ez2, 1/3*ez2, 1/4*phi_2 + 11/12*phi_3 + 3/2*phi_4 + phi_5 - 35/192*ez2; ...
      one, [], [], []; ...
       [], [], [], []; ...
       [], [], [], [] };

c = [ 0; 1/2; 1/2; 1 ];