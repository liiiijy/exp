function s = makestruct(varargin)
% MAKESTRUCT -- Uniformly create input parameter STRUCT.
%               This is merely a utility function to fascilitate
%               easy input syntax to problem-functions.
%
% SYNOPSIS:
%   s = makestruct(struct);
%   s = makestruct('f1', v1, 'f2', v2, ..., 'fn', vn);
%
% PARAMETERS:
%   struct        - Pre-built input parameter struct, returned unchanged
%                   (ie. MAKESTRUCT is nop in this case).
%   'f1', v1, ... - Structure fields 'f1' through 'fn' and corresponding
%                   values v1 through vn.
%
% RETURNS:
%   s - Parmeter structure.
%       Detailed return values as follows:
%       makestruct(struct) -> struct
%       makestruct('f1', v1, ..., 'fn', vn) -> struct('f1', ..., vn)

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.4 $ $Date: 2005/10/12 16:28:00 $

len = numel(varargin);

if len == 1 && isa(varargin{1}, 'struct'),
   s = varargin{1};
elseif len > 1 && mod(len, 2) == 0 && isa(varargin{1}, 'char'),
   s = struct(varargin{:});
else
   error('makestruct:inputwrong', 'Nonsensical input');
end
