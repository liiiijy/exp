function UF = allencahn_post(U, problem)
% ALLENCAHN_POST 
% Postprocessing of numerical solution to Allen--Cahn eqn.
%
% SYNOPSIS:
%   uf = allencahn_post(U, problem);
%
% PARAMETERS:
%   U        - Numerical solution at some (usually final) time step.
%   problem  - Problem specific parameters of the Allen--Cahn problem.
%              Defined by ALLENCAHN. 
%
% RETURNS:
%   UF  - Processed numerical solution `u'.
%
% SEE ALSO:
%   ALLENCAHN, EXPGLM.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.3 $  $Date: 2005/10/22 02:50:13 $


% Ensure x is a row-vector.
x = reshape(problem.x, [1 numel(problem.x)]);

UF = U + x;
