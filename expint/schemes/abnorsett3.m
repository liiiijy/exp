function [u, v, a, b, c] = abnorsett3(z, problem)
% ABNORSETT3 - Adams-Bashforth-Nørsett scheme of stiff order 3.
% 
% SYNOPSIS:
%   [u, v, a, b, c] = abnorsett3(z, problem);
%
% PARAMETERS:
%   z       - Evaluation point (h * linear operator)
%   problem - Problem dependent parameter structure.
%
% RETURNS:
%   u   - The u_ij(z) coefficient functions as a 2-by-3 cell array.
%   a   - The a_ij(z) coefficient functions as a 2-by-2 cell array.
%   v   - The v_ij(z) coefficient function  as a 3-by-3 cell array.
%   b   - The b_ij(z) coefficient functions as a 3-by-2 cell array.
%   c   - The quadrature nodes as a 2-by-1 DOUBLE array.

% This file is part of the 'Expint'-package, 
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.7 $ $Date: 2005/05/11 11:00:14 $

[one, ez2, ez, z] = oneez2(z);

[phi_1, phi_2, phi_3] = phipade(z, 3);

u = { one,                  [],                 []; ...
       ez,  -2*phi_2 - 2*phi_3,  1/2*phi_2 + phi_3 };

a = {                      [],  []; ... 
      phi_1 + 3/2*phi_2+phi_3,  [] };

v = { ez,  -2*phi_2 - 2*phi_3,  1/2*phi_2 + phi_3; ...
      [],                  [],                 []; ...
      [],                 one,                 [] }; 

b = { phi_1 + 3/2*phi_2 + phi_3,  []; ...
                            one,  []; ...
			     [],  [] };

c = [ 0; 1 ];
