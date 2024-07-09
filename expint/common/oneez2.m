function [one, ez2, ez, z] = oneez2(z)
% ONEEZ2 - Helper function for all schemes
%          Determines dimension of the system by checking the
%          input z and returns appropriate numbers.
%
% SYNOPSIS:
%    [one, ez2, ez, z] = oneez2(z)
%
% PARAMETER:
%   z - scalar, vector or matrix which is supposed
%       to represent the the stepsize h times
%       the linear operator of some equation
%
% RETURNS:
%   one - a scalar with one, or a sparse diagonal identity matrix
%   ez2 - exp(z/2) or expm(z/2)
%   ez  - exp(z) or expm(z)
%   z   - the input z, but as a spdiag-matrix if input is vector.
%
% NOTE:
%   For efficiency reasons, ONEEZ2 caches recently computed function
%   values.  The caching behaviour is contingent on the WANTCACHE
%   function and may be toggled on or off as needed.
%
% BUGS:
%   ONEEZ does not support recursive calls on the output z argument.
%   Such calls will treat the incoming z as a matrix, and generate a
%   full version of it.  Protection against this rare case is not
%   implemented as such protection incurs runtime overhead in all cases.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.7 $  $Date: 2005/10/10 08:23:55 $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local user configuration.
atol = 1.0e-10;
rtol = 5.0e-11;

% Uncomment to remove WANTCACHE function dependency (ie. make ONEEZ2 run
% in a self contained environment).
% wantcache = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No user servicable parts below this line.


% Only pay persistent memory overhead if actually needed
if wantcache, persistent cache; end


% Compute new exponentials if the caller forces computation, no previous
% data exists or the existing data is incorrect.

if ~wantcache || isempty(cache) || ~floatequals(cache.z, z, atol, rtol),
   [m, n] = size(z);
   if m == n,
      if m == 1,   %% Scalar case. Simple.
         cache.ez2 = exp(z / 2);
         cache.one = 1;
         cache.z = z;
      else
         %% Matrix case, (square matrix)
         cache.ez2 = expm(full(z) ./ 2);
         cache.one = speye(m, m);
         cache.z = z;
      end
   else
      % Vector case, treated as the main diagonal of a diagonal matrix.
      % Return sparse diagonal matrix to ease programming elsewhere.
      cache.ez2 = spdiags(exp(z ./ 2), 0, m, m);
      cache.one = speye(m, m);
      cache.z   = spdiags(z, 0, m, m);
   end
   cache.ez = cache.ez2 * cache.ez2;
end

% When we get here all we need is in `cache', whether the data was
% (re)computed or not.  Thus, return to caller.

one = cache.one;
ez  = cache.ez;
ez2 = cache.ez2;
z   = cache.z;
