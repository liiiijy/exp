function [u, v, a, b, c] = cfree4(z, problem)
% CFREE4 - Coefficient function for the Commutator-free Lie group 
%          scheme of stiff order 2 with an affine Lie group action.
%
% SYNOPSIS:
%   [u, v, a, b, c] = cfree4(z, problem);
%
% PARAMETERS:
%   z       - Evaluation point (h * linear operator)
%   problem - Problem dependent parameter structure.
%
% RETURNS:
%   u   - The u_ij(z) coefficient functions as a 4-by-1 cell array.
%   a   - The a_ij(z) coefficient functions as a 4-by-4 cell array.
%   v   - The v_ij(z) coefficient function  as a 1-by-1 cell array.
%   b   - The b_ij(z) coefficient functions as a 1-by-4 cell array.
%   c   - The quadrature nodes as a 4-by-1 DOUBLE array. 
%

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.9 $  $Date: 2005/05/11 11:00:14 $

[one, ez2, ez, z] = oneez2(z);

[phi_1]  = phipade(    z, 1);
[phi_12] = phipade(1/2*z, 1);

u = { one; ez2; ez2; ez };

a = {                  [],          [],      [],  []; ...
               1/2*phi_12,          [],      [],  []; ...
                       [],  1/2*phi_12,      [],  []; ...
     1/2*phi_12*(ez2-one),          [],  phi_12,  [] };

v = { ez };

b = { 1/2*phi_1 - 1/3*phi_12,  1/3*phi_1, ...
       1/3*phi_1, -1/6*phi_1 + 1/3*phi_12 };

c = [ 0; 1/2; 1/2; 1 ];
