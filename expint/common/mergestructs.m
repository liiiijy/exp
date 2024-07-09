function s = mergestructs(dflt, ovrd)
% MERGE_STRUCTS -- Merge structures on common and separate fields.
%                  This is a utility function that fascilitate easy
%                  handling of settings supplied to problem-functions.
%
% SYNOPSIS:
%   s = mergestructs(default, override);
%
% PARAMETERS:
%   default  -
%         STRUCT of default values.  Must be non-empty.
%   override -
%         STRUCT of overriding values.
%         These values take precedence over the default values specified
%         in `default'.
%
%         May be empty in which case only default values are used.
%
%         May also contain values for fields unknown to `default'.
%         These values will be set, but only the caller can determine if
%         (and how) such values will be used.
%
% RETURNS:
%   s   - Merged structure with common fields determined by `override'
%         and other fields specified by `default'.
%
% NOTE:
%   This function relies on the `dynamic field names' feature introduced
%   in R13.  Thus, the function does not work in prior releases.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.5 $ $Date: 2005/10/13 13:01:00 $

fn = fieldnames(dflt);
if isempty(fn), 
   error('mergestructs:emptydefault', 'Empty default values...'); 
end

% Set default values.
for idx = 1:numel(fn),
   fld = fn{idx};
   s.(fld) = dflt.(fld);
end

% Reset overriding values.
fn = fieldnames(ovrd);
if ~isempty(fn),
    for idx = 1:numel(fn),
       fld = fn{idx};
       s.(fld) = ovrd.(fld);
    end
end
