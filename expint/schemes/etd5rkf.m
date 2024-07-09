function [u, v, a, b, c] = etd5rkf(z, problem)
% ETD5RKF
% Exponential integrator version based on the fifth order scheme of
% Fehlberg. Six stages.
%
% SYNOPSIS:
%   [u, v, a, b, c] = etd5rkf(z, problem);
%
% PARAMETERS:
%   z       - Evaluation point (h * linear operator)
%   problem - Problem dependent parameter structure.
%
% RETURNS:
%   u   - The u_ij(z) coefficient functions as a 6-by-1 cell array.
%   v   - The b_ij(z)   coefficient function  as a 1-by-1 cell array.
%   a   - The a_ij(z) coefficient functions as a 6-by-6 cell array.
%   b   - The b_ij(z)   coefficient functions as a 1-by-6 cell array.
%   c   - The quadrature nodes as a 6-by-1 DOUBLE array.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.7 $  $Date: 2005/05/11 11:00:14 $

[phi1, phi2, phi3] = phipade(z, 3);
[phi135, phi235]   = phipade(3.0/5.0 .* z, 2);


%%% oneez2 is too specialized for this c-vector, so we copy out code
%%% from there for exponential of our c-vector:
[m, n] = size(z);
if m == n,
   if m == 1,   %% Scalar case. Simple.
      ez29 = exp(2 * z / 9);
      ez13 = exp(    z / 3);
      ez34 = exp(3 * z / 5);
      ez56 = exp(5 * z / 6);
      one = 1;
      ez = exp(z);
   else
      %% Matrix case, (square matrix)
      ez29 = sparse(expm(full(2.*z) ./ 9));
      ez13 = sparse(expm(full(z) ./ 3));
      ez34 = sparse(expm(full(3.*z) ./ 4));
      ez56 = sparse(expm(full(5.*z) ./ 6));
      one  = speye(m, m);
      ez   = sparse(expm(full(z)));
   end
else
   %% Vector case.
   ez29 = spdiags(exp(2 .* z ./ 9), 0, m, m);
   ez13 = spdiags(exp(     z ./ 3), 0, m, m);
   ez34 = spdiags(exp(3 .* z ./ 4), 0, m, m);
   ez56 = spdiags(exp(5 .* z ./ 6), 0, m, m);
   one = speye(m, m);
   ez = spdiags(exp(z), 0, m, m);
end

u = {one; ez29; ez13; ez34; ez; ez56};
c = [0; 2/9; 1/3; 3/4; 1; 5/6];


a = { [], [], [], [], [], []; ...
       -2/3 * phi2 + 10/9 * phi235, [],[],[],[], []; ...
       569/11544 * phi2 + 1355/11544 * phi235, -831/3848 * phi2 + ...
       2755/3848 * phi235, [], [], [], []; ...
   -77157/61568 * phi2 + 143535/61568 * phi235, 587979/61568 * phi2 - ...
       821745/61568 * phi235, ...
       -405/64 * phi2 + 675/64 * phi235, [], [], []; ...
   655263/7696 * phi2 - 2031205/23088 * phi235, -1148769/7696 * phi2 ...
       + 1252665/7696 * phi235, 1593/40 * phi2 - 405/8 * phi235, ...
       144/5 * phi2 - 80/3 * phi235, [], []; ...
   -2212835/277056 * phi2 + 6888625/831168 * phi235, 477285/30784 * phi2 - ...
       496525/30784 * phi235, -39/16 * phi2 + 65/16 ...
       * phi235, -4/9 * phi2 + 20/27 * phi235, -185/96 * phi2 + ...
       575/288 * phi235, [] };

v = { ez };

b = { 47/150 * phi1 - 188/75*phi2 + 94/15 * phi3, ...
       0, ...
       -43/25*phi1 + 132/5*phi2 - 66*phi3, ...
       4124/75 * phi1 - 6152/15 * phi2 + 2704/3 * phi3, ...
       189/10 * phi1 - 662/5 * phi2 + 284 * phi3, ...
       -1787/25 * phi1 + 12966/25 * phi2 - 5628/5 * phi3 };
