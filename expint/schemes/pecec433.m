function [u, v, a, b, c] = pecec433(z, problem)
% PECEC433 - Order 4 exponential general linear method 
%            with s=3 stages and r=3 outputs.
% 
% SYNOPSIS:
%   [u, v, a, b, c] = pecec433(z, problem);
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

[phi_1, phi_2, phi_3, phi_4] = phipade(z, 4, 13);

u = { one, [], []; ...
      ez, -2*phi_2 - 2*phi_3, 1/2*phi_2 + phi_3; ...
      ez, -phi_2 + phi_3 + 3*phi_4, 1/6*phi_2 - phi_4 };

a = {                                    [],                        [], []; ...
                  phi_1 + 3/2*phi_2 + phi_3,                        [], []; ...
      phi_1 + 1/2*phi_2 - 2*phi_3 - 3*phi_4, 1/3*phi_2 + phi_3 + phi_4, [] };

v = { ez, -phi_2 + phi_3 + 3*phi_4, 1/6*phi_2 - phi_4; ...
      [],  [], []; ...
      [], one, [] };

b = { phi_1 + 1/2*phi_2 - 2*phi_3 - 3*phi_4, [], 1/3*phi_2 + phi_3 + phi_4; ...
      one, [], []; ...                      
       [], [], [] };

c = [ 0; 1; 1 ];
