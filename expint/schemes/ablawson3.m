function [u, v, a, b, c] = ablawson3(z, problem)
% ABLAWSON3 - Adams-Bashforth-Lawson scheme of stiff order 1. 
%
% SYNOPSIS:
%   [u, v, a, b, c] = ablawson3(z, problem);
%
% PARMETERS:
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
% $Revision: 1.5 $  $Date: 2005/10/10 07:22:11 $

[one, ez2, ez, z] = oneez2(z);

e2z = ez*ez;
e3z = e2z*ez;

u = { one,        [],        []; ...
       ez,  -4/3*e3z,  5/12*e3z };

a = {        [],  []; ...
       23/12*ez,  [] }; 

v = { ez, -4/3*e2z, 5/12*e3z; ...
      [],       [],       []; ... 
      [],      one,       []};

b = { 23/12*ez, []; ...
           one, []; ...
            [], [] };
    
c = [ 0; 1 ];
