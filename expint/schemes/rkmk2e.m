function [u, v, a, b, c] = rkmk2e(z, problem)
% RKMK2E - Coefficient function for the Runge-Kutta-Munthe-Kaas scheme
%          with affine Lie group action based on the midpoint rule, it
%          has stiff order 2.
%
% SYNOPSIS:
%   [u, v, a, b, c] = rkmk2e(z, problem);
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

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.6 $  $Date: 2005/10/10 07:22:12 $

[one, ez2, ez, z] = oneez2(z);

phi_1 = phipade(z, 1);

u = { one;  ez };

a = {    [],    []; 
      phi_1,    [] };

v = { ez };

b = { 1/2*phi_1,  1/2*phi_1 };

c = [ 0; 1 ];
