function [u, v, a, b, c] = lawson2b(z, problem)
% LAWSON2B - Runge--Kutta--Lawson based on the trapeziodal rule with
%            stiff order 1.
%
% SYNOPSIS:
%   [u, v, a, b, c] = lawson2b(z, problem);
%
% PARAMETERS:
%   z       - Evaluation point (h * linear operator)
%   problem - Problem dependent parameter structure.
%
% RETURNS:
%   u   - The u_ij(z) coefficient functions as a 2-by-1 cell array.
%   a   - The a_ij(z) coefficient functions as a 2-by-2 cell array.
%   v   - The v_ij(z) coefficient function  as a 1-by-1 cell array.
%   b   - The b_ij(z) coefficient functions as a 1-by-2 cell array.
%   c   - The quadrature nodes as a 2-by-1 DOUBLE array.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.7 $  $Date: 2005/10/10 07:22:12 $

[one, ez2, ez, z] = oneez2(z);

u = { one; ez };

a = { [], []; ... 
      ez, [] };

v = { ez };

b = { 1/2*ez, 1/2*one };

c = [ 0; 1 ];
