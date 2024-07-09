function [u, v, a, b, c] = abnorsett4(z, problem)
% ABNORSETT4 - Adams-Bashforth-Nørsett scheme of stiff order 4.
% 
% SYNOPSIS:
%   [u, v, a, b, c] = abnorsett4(z, problem);
%
% PARAMETERS:
%   z       - Evaluation point (h * linear operator)
%   problem - Problem dependent parameter structure.
%
% RETURNS:
%   u   - The u_ij(z) coefficient functions as a 2-by-4 cell array.
%   a   - The a_ij(z) coefficient functions as a 2-by-2 cell array.
%   v   - The v_ij(z) coefficient function  as a 4-by-4 cell array.
%   b   - The b_ij(z) coefficient functions as a 4-by-2 cell array.
%   c   - The quadrature nodes as a 2-by-1 DOUBLE array.

% This file is part of the 'Expint'-package, 
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.11 $ $Date: 2005/06/21 17:37:48 $

[one, ez2, ez, z] = oneez2(z);

[phi_1, phi_2, phi_3, phi_4] = phipade(z, 4, 13);

u = {one,    [],    [],    []; ...
     ez,   -3*phi_2 - 5*phi_3 - 3*phi_4, ...
     3/2*phi_2 + 4*phi_3 + 3*phi_4,   -1/3*phi_2 - phi_3 - phi_4 };

a = {                                    [],  []; ...
       phi_1 + 11/6*phi_2 + 2*phi_3 + phi_4,  [] };

v = { ez, -3*phi_2 - 5*phi_3 - 3*phi_4, ...
      3/2*phi_2 + 4*phi_3 + 3*phi_4, -1/3*phi_2 - phi_3 - phi_4; ...
      [],  [],  [], []; ...
      [], one,  [], []; ...
      [],  [], one, [] }; 

b = {phi_1 + 11/6*phi_2 + 2*phi_3 + phi_4, []; ...
                                      one, []; ...
                                       [], []; ...
                                       [], [] };

c = [ 0; 1 ];
