function [u, v, a, b, c] = eglm433(z, problem)
% EGLM433 - Order 4 exponential general linear method 
%           with s=3 stages and r=3 outputs.
% 
% SYNOPSIS:
%   [u, v, a, b, c] = eglm433(z, problem);
%
% PARAMETERS:
%   z       - Evaluation point (h * linear operator)
%   problem - Problem dependent parameter structure.
%
% RETURNS:
%   u   - The u_ij(z) coefficient functions as a 3-by-3 cell array.
%   a   - The a_ij(z) coefficient functions as a 3-by-3 cell array.
%   v   - The v_ij(z) coefficient function  as a 3-by-3 cell array.
%   b   - The b_ij(z) coefficient functions as a 3-by-3 cell array.
%   c   - The quadrature nodes as a 3-by-1 DOUBLE array.

% This file is part of the 'Expint'-package, 
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.1 $ $Date: 2005/06/21 17:37:48 $

[one, ez2, ez, z] = oneez2(z);

[phi_1, phi_2, phi_3, phi_4, phi_5] = phipade(z, 5, 13);
[phi_12, phi_22, phi_32] = phipade(1/2*z, 3, 13);

u = { one, [], []; ...
      ez2, -1/2*phi_22-1/4*phi_32, 1/8*phi_22 + 1/8*phi_32; ...
      ez, -2/3*phi_2 + 2*phi_3 + 4*phi_4, 1/10*phi_2 - 1/5*phi_3 - 6/5*phi_4 };

a = {                      [],               [], []; ...
       1/2*phi_12 + 3/8*phi_22 + 1/8*phi_32, [], []; ...
      phi_1 - 1/2*phi_2 - 5*phi_3 - 6*phi_4, ...
        16/15*phi_2 + 16/5*phi_3 + 16/5*phi_4, [] };

v = { ez, -1/3*phi_2 + 5/3*phi_3 - phi_4 - 8*phi_5, ...
      1/30*phi_2 - 2/15*phi_3 - 1/5*phi_4 + 8/5*phi_5; ...
      [],  [], []; ...
      [], one, [] };

b = { phi_1 - 3/2*phi_2 - 4*phi_3 + 9*phi_4 + 24*phi_5, ...
      32/15*phi_2 + 32/15*phi_3 - 64/5*phi_4 - 128/5*phi_5, ...
      -1/3*phi_2 + 1/3*phi_3 + 5*phi_4 + 8*phi_5; ...
                                     one, [], []; ...
                                      [], [], [] };

c = [ 0; 1/2; 1 ];
