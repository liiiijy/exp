function [u, v, a, b, c] = abnorsett2(z, problem)
% ABNORSETT2 - Adams-Bashforth-Nørsett scheme of stiff order 2.
% 
% SYNOPSIS:
%   [u, v, a, b, c] = abnorsett2(z, problem);
%
% PARAMETERS:
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
% $Revision: 1.6 $ $Date: 2005/05/11 08:30:27 $

[one, ez2, ez, z] = oneez2(z);

[phi_1, phi_2] = phipade(z, 2);

u = { one,     []; ...
       ez, -phi_2 };

a = {             [],  []; ...
       phi_1 + phi_2,  [] };

v = {  ez, -phi_2; ...
       [],     [] };

b = { phi_1 + phi_2, []; ...
                one, [] };

c = [ 0; 1 ];
