function [u, v, a, b, c] = etd3rk(z, problem)
% ETD3RK - The ETD3RK coefficent functions
%     
% SYNOPSIS:
%   [u, v, a, b, c] = etd3rk(z, problem);
%
% PARAMETERS:
%   z       - Evaluation point (h * linear operator)
%   problem - Problem dependent parameter structure.
%
% RETURNS:
%   u   - The u_ij(z) coefficient functions as a 3-by-1 cell array.
%   a   - The a_ij(z) coefficient functions as a 3-by-3 cell array.
%   v   - The v_ij(z) coefficient function  as a 1-by-1 cell array.
%   b   - The b_ij(z) coefficient functions as a 1-by-3 cell array.
%   c   - The quadrature nodes as a 3-by-1 DOUBLE array.


% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.8 $  $Date: 2005/05/11 16:17:00 $

[one, ez2, ez, z] = oneez2(z);

[phi_1, phi_2, phi_3] = phipade(z, 3);
[phi_12]              = phipade(1/2*z, 1);

u = { one; ez2; ez };

a  = {          [],          [],     []; ...
        1/2*phi_12,          [],     []; ...
            -phi_1,     2*phi_1,     [] };


v = { ez };

b = {   phi_1 - 3*phi_2 + 4*phi_3, ...
       4*phi_2 - 8*phi_3, ...
       -phi_2 + 4*phi_3 };

c = [ 0; 1/2; 1 ];
