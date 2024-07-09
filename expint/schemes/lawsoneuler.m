function [u, v, a, b, c] = lawsoneuler(z, problem)
% LAWSONEULER - Lawson--Euler scheme of stiff order 1.
%
% SYNOPSIS:
%   [u, v, a, b, c] = lawsoneuler(z, problem);
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
% $Revision: 1.5 $  $Date: 2005/10/10 07:22:12 $

[one, ez2, ez, z] = oneez2(z);

u = { one };

a = { [] };

v = { ez };

b = { ez };

c = [ 0 ];
