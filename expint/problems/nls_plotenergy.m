function nls_plotenergy(u, t,  problem)
% NLS_PLOTENERGY -- Output plotting function for the NLS problem.
%
% SYNOPSIS:
%   nls_plotenergy(u, t, problem)
%
% DESCRIPTION:
%   Plots the energy of the solution as a function of time.
%
% PARAMETERS:
%   u        - (Numerical) solution to the NLS problem.
%   t        - point in time. 
%   problem  - Problem description structure as defined by NLS.
%
% RETURNS:
%   Nothing.
%
% SEE ALSO:
%   NLS, EXPGLM.
%
% LIMITATIONS:
%      - Uses global variables.


% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.2 $  $Date: 2005/10/22 02:50:13 $

% Variables we need here must also be declared global
% in the calling environment for this to work!
% You should also remember to reset all variables between each run.

global energydata timeindex timedata
% (again: repeat this line in your calling environment)

if (sum(timeindex) == 0 | isempty(timeindex)),
   timeindex = 1; 
end

timedata(timeindex) = t;

energydata(timeindex) = nls_energy(u, pr);

timeindex = timeindex + 1;

plot(timedata, energydata);

drawnow


