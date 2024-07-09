function [u, v, a, b, c] = modgenlawson42(z, problem)
% MODGENLAWSON42 - Modified generalized Lawson scheme of stiff order 3.
%
% SYNOPSIS:
%   [u, v, a, b, c] = modgenlawson42(z, problem);
%
% PARAMETERS:
%   z       - Evaluation point (h * linear operator)
%   problem - Problem dependent parameter structure.
%
% RETURNS:
%   u   - The u_ij(z) coefficient functions as a 4-by-2 cell array.
%   a   - The a_ij(z) coefficient functions as a 4-by-4 cell array.
%   v   - The v_ij(z) coefficient function  as a 2-by-2 cell array.
%   b   - The b_ij(z) coefficient functions as a 2-by-4 cell array.
%   c   - The quadrature nodes as a 4-by-1 DOUBLE array.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.8 $ $Date: 2005/05/11 11:00:14 $

[one, ez2, ez, z] = oneez2(z);

[phi_1, phi_2, phi_3] = phipade(z, 3);
[phi_12, phi_22]      = phipade(1/2*z, 2);

u = { one,                    []; ...
      ez2,           -1/4*phi_22; ...
      ez2, -1/4*phi_22 + 1/4*one; ...
       ez,      -phi_2 + 1/2*ez2 };
 
a = {                               [],      [],   [],  []; ...
               1/2*phi_12 + 1/4*phi_22,      [],   [],  []; ...
     1/2*phi_12 + 1/4*phi_22 - 3/4*one, 1/2*one,   [],  []; ...
               phi_1 + phi_2 - 3/2*ez2,      [],  ez2,  [] };
 
v = { ez, -1/2*phi_2 + phi_3 + 1/12*ez2; ...
      [],                            [] };

b = { phi_1 - 2*phi_3 - 1/2*ez2, 1/3*ez2, 1/3*ez2, 1/2*phi_2 + phi_3 - 1/4*ez2; ...
                            one,      [],      [],                [] };

c = [ 0; 1/2; 1/2; 1 ];
