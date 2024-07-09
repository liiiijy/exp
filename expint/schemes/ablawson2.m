function [u, v, a, b, c] = ablawson2(z, problem)
% ABLAWSON2 - Adams-Bashforth-Lawson scheme of stiff order 1.
% 
% SYNOPSIS:
%   [u, v, a, b, c] = ablawson2(z, problem);
%
% PARMETERS:
%   z       - Evaluation point (h * linear operator)
%   problem - Problem dependent parameter structure.
%
% RETURNS:
%   u   - The u_ij(z) coefficient functions as a 2-by-2 cell array.
%   a   - The a_ij(z) coefficient functions as a 2-by-2 cell array.
%   v   - The v_ij(z) coefficient function  as a 2-by-2 cell array.
%   b   - The b_ij(z) coefficient functions as a 2-by-2 cell array.
%   c   - The quadrature nodes as a 2-by-1 DOUBLE array.

% This file is part of the 'Expint'-package, 
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.4 $  $Date: 2005/05/11 08:30:27 $

[one, ez2, ez, z] = oneez2(z);

e2z = ez*ez;

u = {  one,       []; ... 
        ez,  -1/2*e2z };

a = {      [],  []; ...
       3/2*ez,  [] };

v = {  ez, -1/2*ez2; ...
       [],       [] }; 

b = { 3/2*ez, []; ... 
         one, [] };

c = [ 0; 1 ];
