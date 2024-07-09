function [u, v, a, b, c] = norsetteuler(z, problem)
% NORSETTEULER - Nørsett-Euler scheme of stiff order 1.
%
%   Alias names: ETD-Euler, Lie-Euler, ETD1RK, exponentially
%                fitted Euler, filtered Euler.
%
% SYNOPSIS:
%   [u, v, a, b, c] = norsetteuler(z, problem);
%
% PARAMETERS:
%   z       - Evaluation point (h * linear operator)
%   problem - Problem dependent parameter structure.
%
% RETURNS:
%   u   - The u_ij(z) coefficient functions as a 1-by-1 cell array.
%   a   - The a_ij(z) coefficient functions as a 1-by-1 cell array.
%   v   - The v_ij(z) coefficient function  as a 1-by-1 cell array.
%   b   - The b_ij(z) coefficient functions as a 1-by-1 cell array.
%   c   - The quadrature nodes as a 1-by-1 DOUBLE array.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.2 $  $Date: 2005/10/10 07:22:12 $

[one, ez2, ez, z] = oneez2(z);

phi_1 = phipade(z, 1);

u = { one };

a = { [] };

v = { ez };

b = { phi_1 };

c = [ 0 ];
