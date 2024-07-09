function [u, v, a, b, c] = pec524(z, problem)
% PEC524 - Order 5 exponential general linear method 
%           with s=2 stages and r=4 outputs.
% 
% SYNOPSIS:
%   [u, v, a, b, c] = pec524(z, problem);
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
% $Revision: 1.1 $ $Date: 2005/06/24 12:31:49 $

[one, ez2, ez, z] = oneez2(z);

[phi_1, phi_2, phi_3, phi_4, phi_5] = phipade(z, 5, 13);

u = { one, [], [], []; ...
      ez, -3*phi_2 - 5*phi_3 - 3*phi_4, ...
      3/2*phi_2 + 4*phi_3 + 3*phi_4, ...
      -1/3*phi_2 - phi_3 - phi_4 };

a = {                                    [],  []; ...
       phi_1 + 11/6*phi_2 + 2*phi_3 + phi_4,  [] };

v = { ez, -3/2*phi_2 + 1/2*phi_3 + 6*phi_4 + 6*phi_5, ...
      1/2*phi_2 + 1/3*phi_3 - 3*phi_4 - 4*phi_5, ...
      -1/12*phi_2 - 1/12*phi_3 + 1/2*phi_4 + phi_5; ...
      [],  [],  [],  []; ...
      [], one,  [],  []; ...
      [],  [], one,  [] }; 

b = { phi_1 + 5/6*phi_2 - 5/3*phi_3 - 5*phi_4 - 4*phi_5, ...
      1/4*phi_2 + 11/12*phi_3 + 3/2*phi_4 + phi_5; ...
                                    one, []; ...
                                     [], []; ...
                                     [], [] };

c = [ 0; 1 ];
