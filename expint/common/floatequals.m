function B = floatequals(x, y, varargin)
% FLOATEQUALS - Determine floating point equality for arrays of DOUBLE or
%          COMPLEX.
%
% SYNOPSIS:
%   B = floatequals(x, y);
%   B = floatequals(x, y, atol);
%   B = floatequals(x, y, atol, rtol);
%
% DESCRIPTION:
%   Using direct equality operators (== or ~=) may not be appropriate for
%   matrices of class DOUBLE or COMPLEX.  FLOATEQUALS implements a weaker
%   sense of equality with exact equality definitions overridable by the
%   user.
%
%   Two objects, X and Y, are deemed equal if and only if
%    - ALL(SIZE(X) == SIZE(Y)), and
%    - (ALL(ABS(X - Y) < atol) or
%       ALL(ABS(X - Y) < rtol.*ABS(Y)))
%   with `atol' and `rtol' being absolute and relative tolerances
%   respectively.
%
%   For complex arrays, separate checks are made for the real and
%   imaginary parts.
%
% PARAMETERS:
%   x, y - Objects to check for equality.
%   atol - Absolute tolerance.  OPTIONAL.  DEFAULT VALUE = 1.0e-6.
%   rtol - Relative tolerance.  OPTIONAL.  DEFAULT VALUE = 1.0e-7.
%
% RETURNS:
%   B    - Boolean status indicating whether x is equal to y or not.
%          Possible values are  TRUE  and  FALSE.
%
% SEE ALSO:
%   RELOP, REAL, IMAG, TRUE, FALSE.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.5 $ $Date: 2005/05/13 18:07:48 $

% DUPLICATION NOTE:
%
%   This function is duplicated as FLTEQ in PHIPADE in order to make
%   PHIPADE self contained.  Any change made to this function should be
%   repeated in PHIPADE:FLTEQ.

error(nargchk(2, 4, nargin));

if issparse(x) || issparse(y),
   % only work on non-zero elements...
   [ix, jx, x] = find(x);
   [iy, jy, y] = find(y);

   % x ~= y if different sparsity patterns
   B = (numel(ix) == numel(iy)) && ...
       all(ix == iy) && all(jx == jy);
else
   B = true;                    % additional testing needed
end

sx = size(x);
sy = size(y);

if B && all(sx == sy),
   [atol, rtol] = deal(1.0e-6, 1.0e-7);

   if nargin > 2, atol = abs(varargin{1}); end
   if nargin > 3, rtol = abs(varargin{2}); end

   % Straighten out  x  and  y  for multi-D cases
   if all(sx > 1), x = x(:); y = y(:); end

   [xc, yc] = deal(real(x), real(y));
   xc = abs(xc - yc);
   a = all(xc < atol);
   r = all(xc < rtol.*abs(yc));

   if ~isreal(x) || ~isreal(y),
      [xc, yc] = deal(imag(x), imag(y));
      xc = abs(xc - yc);
      a = a & all(xc < atol);
      r = r & all(xc < rtol.*abs(yc));
   end

   B = a | r;
else
   B = false;
end
