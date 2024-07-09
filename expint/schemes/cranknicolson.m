function [u, v, a, b, c] = cranknicolson(z, problem)
% CRANKNICOLSON - Coefficient functions for an implementation of
%                 Crank-Nicolson in a exponential integrator framework.
%
%   The scheme is implicit, and thus a Newton iterator is
%   implemented. Only the diagonal part is included in the inverse which
%   is formed to make it computationally realistic. 
%
%   Crank--Nicolson is the method
%   
%      y1 = y0 + h/2 * (f(y0) + f(y1)) where f(y0) = A*y0 + b(y0)
%
%   The equation to solve by Newtons method becomes
%   
%      F(y1) = 0   where
%         F(y1) = y1 - h/2 A*y1 - h/2 b(y1) - yn + h/2 A*yn + h/2 b(yn)
%      where y1 is our candidate solution and yn is some iterate (equal
%      to y0 in the first iteration).
%
%   The Newton iteration is then 
%  
%      y1 = yn - F'^(-1)(yn) F(yn)  where the b'(yn) part of F' is removed
%                                   when inverting
%
%   which becomes after inserting the more simple
%  
%      y1 = F'^(-1) * (h/2 b(yn) + y0 + h/2 A y0 + h/2 b(y0))
%
%   Limiting ourselve to four iterations, the above fits nicely
%   into the format of an exponential integrator
%


% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.5 $  $Date: 2005/05/11 11:00:14 $


% We only use 'one' from this function, but usually it is not
% an overhead since the data is usually persistent.
[one, ez2, ez, z] = oneez2(z);

pade11 = (one - z./2) \ (one + z./2);
phiCN = (2*one - z) \ one;

u = {one; pade11;  pade11; pade11};

a = {       [],     [],     [],  []; ...
      2.*phiCN,     [],     [],  []; ...
         phiCN,  phiCN,     [],  []; ...
         phiCN,     [],  phiCN,  [] };

v = { pade11 };

b = { phiCN,  [],  [],  phiCN };

c = [ 0; 1; 1; 1 ];
