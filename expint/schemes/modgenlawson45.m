function [u, v, a, b, c] = modgenlawson45(z, problem)
% MODGENLAWSON45 - Modified generalized Lawson scheme of stiff order 6.
%
% SYNOPSIS:
%   [u, v, a, b, c] = modgenlawson45(z, problem);
%
% PARAMETERS:
%   z       - Evaluation point (h * linear operator)
%   problem - Problem dependent parameter structure.
%
% RETURNS:
%   u   - The u_ij(z) coefficient functions as a 4-by-5 cell array.
%   a   - The a_ij(z) coefficient functions as a 4-by-4 cell array.
%   v   - The v_ij(z) coefficient function  as a 5-by-5 cell array.
%   b   - The b_ij(z) coefficient functions as a 5-by-4 cell array.
%   c   - The quadrature nodes as a 4-by-1 DOUBLE array.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.8 $ $Date: 2005/05/19 17:14:03 $

[one, ez2, ez, z] = oneez2(z);

[phi_1, phi_2, phi_3, phi_4, phi_5, phi_6]  = phipade(z, 6);
[phi_12, phi_22, phi_32, phi_42, phi_52]    = phipade(1/2*z, 5);

u = { one, [], [], [], []; ...
      ez2, -phi_22 - 13/12*phi_32 - 9/16*phi_42 - 1/8*phi_52, ...
      3/4*phi_22 + 19/16*phi_32 + 3/4*phi_42 + 3/16*phi_52, ...
     -1/3*phi_22 - 7/12*phi_32 - 7/16*phi_42 - 1/8*phi_52, ...
      1/16*phi_22 + 11/96*phi_32 + 3/32*phi_42 + 1/32*phi_52; ...
      ez2, -phi_22 - 13/12*phi_32 - 9/16*phi_42 - 1/8*phi_52 + 105/64*one, ...
      3/4*phi_22 + 19/16*phi_32 + 3/4*phi_42 + 3/16*phi_52 - 189/128*one, ...
     -1/3*phi_22 - 7/12*phi_32 - 7/16*phi_42 - 1/8*phi_52 + 45/64*one, ...
      1/16*phi_22 + 11/96*phi_32 + 3/32*phi_42 + 1/32*phi_52 - 35/256*one; ...
      ez, -4*phi_2 - 26/3*phi_3 - 9*phi_4 - 4*phi_5 + 105/32*ez2, ...
      3*phi_2 + 19/2*phi_3 + 12*phi_4 + 6*phi_5 - 189/64*ez2, ...
     -4/3*phi_2 - 14/3*phi_3 - 7*phi_4 - 4*phi_5 + 45/32*ez2, ...
      1/4*phi_2 + 11/12*phi_3 + 3/2*phi_4 + phi_5 - 35/128*ez2 };

a = { [], [], [], []; ...
      1/2*phi_12 + 25/48*phi_22 + 35/96*phi_32 + ...
      5/32*phi_42 + 1/32*phi_52, [], [], []; ...
      1/2*phi_12 + 25/48*phi_22 + 35/96*phi_32 + ...
      5/32*phi_42 + 1/32*phi_52 - 315/256*one, 1/2*one, [], []; ...
      phi_1 + 25/12*phi_2 + 35/12*phi_3 + 5/2*phi_4 + ...
      phi_5 - 315/128*ez2, [], ez2, [] };

b = { -3055/3776 * ez2 - 935/708 * phi_3 + 300/59 * phi_6 + ...
      755/708 * phi_2 - 755/118 * phi_4 - 541/59 * phi_5 + phi_1, ...
      1/3 * ez2, 1/3 * ez2, 50/59 * phi_3 + 12/59 * phi_2 - ...
      60/59 * phi_6 - 157/944 * ez2 + 120/59 * phi_5 + 105/59 * phi_4; ...
      one, [], [], []; ...
      [], [], [], []; ...
      [], [], [], []; ...
      [], [], [], [] };

v = { ez, 495/944 * ez2 - 34/177 * phi_3 - 600/59 * phi_6 - ...
      116/59 * phi_2 + 519/59 * phi_4 + 964/59 * phi_5, ...
      -577/1888 * ez2 + 57/59 * phi_2 + 121/118 * phi_3 + ...
      600/59 * phi_6 - 342/59 * phi_4 - 846/59 * phi_5, ... 	
      25/236 * ez2 - 76/177 * phi_3 - 300/59 * phi_6 - ...
      56/177 * phi_2 + 112/59 * phi_4 + 364/59 * phi_5, ...
      -181/11328 * ez2 + 49/708 * phi_3 + 60/59 * phi_6 + ...
      11/236 * phi_2 - 33/118 * phi_4 - 61/59 * phi_5; ...
      [], [], [], [], []; ...
      [], one, [], [], []; ...
      [], [], one, [], []; ...
      [], [], [], one, [] };

c = [ 0; 1/2; 1/2; 1 ];
