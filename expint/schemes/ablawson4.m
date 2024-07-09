function [u, v, a, b, c] = ablawson4(z, problem)
% ABLAWSON4 - Adams-Bashforth-Lawson scheme of stiff order 1.
% 
% SYNOPSIS:
%   [u, v, a, b, c] = ablawson4(z, problem);
%
% PARMETERS:
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
% $Revision: 1.5 $  $Date: 2005/05/11 11:00:14 $

[one, ez2, ez, z] = oneez2(z);

e2z = ez*ez;
e3z = e2z*ez;
e4z = e2z*e2z;
    
u = { one,         [],        [],       []; ...
       ez, -59/24*e2z, 37/24*e3z, -3/8*e4z };

a = {        [], []; ...
       55/24*ez, [] };

v = { ez, -59/24*e2z, 37/24*e3z, -3/8*e4z; ...
      [],         [],        [],       []; ...
      [],        one,        [],       []; ...
      [],         [],       one,       [] };

b = { 55/24*ez, []; ... 
           one, []; ...
	    [], []; ...
	    [], [] };

c = [ 0; 1 ];
