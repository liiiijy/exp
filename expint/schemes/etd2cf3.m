function [u, v, a, b, c] = etd2cf3(z, problem)
% ETD2CF3 - Coefficient for a third order ETD version based on a 
%           commutator-free scheme by Celledoni, Martinsen and Owren.
%     
% SYNOPSIS: 
%  [u, v, a, b, c] = ETD2CF3(z, problem)
%
% PARAMETERS:
%   z       - Evaluation point (h * linear operator)
%   problem - Problem dependent parameter structure.
%
% RETURNS:
%   u   - The u_ij(z) coefficient functions as a 3-by-1 cell array.
%   a   - The a_ij(z) coefficient functions as a 3-by-3 cell array.
%   v   - The v_ij(z) coefficient function  as a 1-by-1 cell array.
%   b   - The b_ij(z) coefficient functions as a 1-by-3 cell array.
%   c   - The quadrature nodes as a 3-by-1 DOUBLE array.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.4 $  $Date: 2005/05/11 16:29:22 $

[one, ez2, ez, z] = oneez2(z);

[phi_1, phi_2, phi_3] = phipade(z, 3);
 phi_12               = phipade(1/2*z, 1);
[phi_13, phi_23]      = phipade(2/3*z, 2);

%%% oneez2 is too specialized for this c-vector, so we copy out code
%%% from there for exponential of our c-vector:
[m, n] = size(z);
if m == n,
   if m == 1,   %% Scalar case. Simple.
      one = 1;
      ez13 = exp(1/3*z);
      ez23 = exp(2/3*z);
      ez = exp(z);
   else
      %% Matrix case, (square matrix)
      one  = speye(m, m);
      ez13 = sparse(expm(1/3*full(z)));
      ez23 = sparse(expm(1/3*full(2.*z)));
      ez   = sparse(expm(full(z)));
   end
else
   %% Vector case.
   one = speye(m, m);
   ez13 = spdiags(exp(1/3*z), 0, m, m);
   ez23 = spdiags(exp(2/3*z), 0, m, m);
   ez = spdiags(exp(z), 0, m, m);
end

u = { one; ez13; ez23 };

a  = {                      [],             [],     []; ...
                    1/3*phi_12,             [],     []; ...
       2/3*phi_13 - 4/3*phi_23,     4/3*phi_23,     [] };

v = { ez };

b = {  phi_1 - 9/2*phi_2 + 9*phi_3, ...
       6*phi_2 - 18*phi_3, ...
       -3/2 * phi_2 + 9*phi_3 };

c = [ 0; 1/3; 2/3 ];
