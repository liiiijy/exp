function M = setupschemes(varargin)
% SETUPSCHEMES - Function to set up a cell-array of
%                schemestructure defining schemes to be
%                used for testing.
%
% SYNOPSIS:
%    M = setupschemes(s1, s2, ..., sn);
%
% PARAMETERS:
%   s<k> - Name (string) of coefficient function <k>.
%          At least one scheme name must be provided.
%
% RETURNS:
%   M - Structure cell-array, one struct for each scheme.
%       The structure holds the following information
%           coffcn - the matlab function for the scheme (string)
%           name - the name of the scheme, equals `coffcn', may be used
%                  in plot legends.
%           linestyle - linestyles picked from a list of linestyles.
%
%   Other scripts may make use of the following optional fields:
%           workscaling - a number that could be used by anyone in order
%                  to scale comparisons among schemes.
%           relstages - Set to a number different from 1 in order to
%                  divide the stepsize used by this number. If you
%                  are comparing some schemes with others, you might
%                  want to do this in order to make the comparison fair.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.13 $  $Date: 2005/10/12 16:28:00 $

error(nargoutchk(0, 2, nargout));

default_linestyles = { 'k.-',  'y+-',  'ro-',  'bx-',  'g*-',  'cs-',  'md-',  ...
                       'k.--', 'y+--', 'ro--', 'bx--', 'g*--', 'cs--', 'md--', ...
                       'k.-.', 'y+-.', 'ro-.', 'bx-.', 'g*-.', 'cs-.', 'md-.' };

if nargin > 0,
   args = varargin;
   narg = 0;

   for i = 1:numel(args),
      if which(args{i}),        % function is in path (ie. callable)
         j = mod(narg, numel(default_linestyles)) + 1;

         M{narg+1} = struct('coffcn', args{i}, ...
                            'name', args{i},   ...
                            'linestyle',       ...
                             default_linestyles{j});
         narg = narg + 1;
      else
         warning('setupscheme:notfound', ...
                 ['Function `', args{i}, ''' not in path']);
      end
   end

   if narg == 0,
      error('setupschemes:invalid', ...
            'No valid stepping schemes supplied');
   end
else
   error('setupschemes:noscheme', 'No stepping scheme supplied');
end
