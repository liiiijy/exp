function [u, v, a, b, c] = eglm332(z, problem)
% EGLM332 - Order 4 exponential general linear method 
%           with s=2 stages and r=3 outputs.
% 
% SYNOPSIS:
%   [u, v, a, b, c] = eglm332(z, problem);
%
% PARAMETERS:
%   z       - Evaluation point (h * linear operator)
%   problem - Problem dependent parameter structure.
%
% RETURNS:
%   u   - The u_ij(z) coefficient functions as a 3-by-2 cell array.
%   a   - The a_ij(z) coefficient functions as a 3-by-3 cell array.
%   v   - The v_ij(z) coefficient function  as a 2-by-2 cell array.
%   b   - The b_ij(z) coefficient functions as a 2-by-3 cell array.
%   c   - The quadrature nodes as a 3-by-1 DOUBLE array.

% This file is part of the 'Expint'-package, 
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.1 $ $Date: 2005/06/21 17:37:48 $

[one, ez2, ez, z] = oneez2(z);

[phi_1, phi_2, phi_3, phi_4] = phipade(z, 4, 13);
[phi_12, phi_22] = phipade(1/2*z, 2, 13);

u = { one, []; ...
      ez2, -1/4*phi_22; ...
      ez, -1/3*phi_2 + 4/3*phi_3 };

a = {                      [],                     [], []; ...
      1/2*phi_12 + 1/4*phi_22,                     [], []; ...
      phi_1 - phi_2 - 4*phi_3,  4/3*phi_2 + 8/3*phi_3, [] };

v = { ez, -1/6*phi_2 + phi_3 - 2*phi_4; ...
      [],  [] };

b = { phi_1 - 2*phi_2 - 2*phi_3 + 12*phi_4, ...
                      8/3*phi_2 - 16*phi_4, ...
              -1/2*phi_2 + phi_3 + 6*phi_4; ...
                               one, [], [] };

c = [ 0; 1/2; 1 ];
