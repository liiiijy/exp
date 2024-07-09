function H = sinegordon_hamilt(y, problem)
% SINEGORDON_HAMILT - Hamiltonian function for the sine-Gordon equation.
%
% SYNOPSIS:
%   H = sinegordon_hamilt(y, problem);
% 
% DESCRIPTION:
%   Evaluates the Hamiltonian function
%     H(p, q) = \int_0^Lengthp^2/2 + q_x^2/2 + 1 - cos(q) dx
%
%   This function may also be specified as the post processor of
%   the sine-Gordon problem.
%
% PARAMETERS:
%   y       - Evaluation point in Fourier space, y = (q, p).
%             If input is matrix, each row is a solution value for
%             which the Hamiltonian is computed.
%   problem - Problem dependent parameters.  Should be defined by SINEGORDON.
%
% RETURNS:
%   H      - The value(s) of the Hamiltonian function at the given
%            phase value(s).
%
% SEE ALSO:
%   SINEGORDON, SINEGORDON_POST

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.1 $  $Date: 2005/11/09 03:18:29 $

ND = problem.ND;

% If not matrix input:
if min(size(y)) == 1
   % Reshape into row-vector.
   y = reshape(y, 1, 2*ND);
end

% From here we assume that y(n,:) (each row) is a solution at time t_n.
numsolutions = size(y, 1);
H = zeros(numsolutions, 1);

dx = problem.Length/ND;

for n=1:numsolutions,
   ys = reshape(y(n,:), 2*ND, 1);
   
   q = ifft(ys(1:ND));
   p = ifft(ys(ND+1:2*ND));
   
   % Calculate spatial derivative, first order difference.
   qx = 1/dx * diff([q ; q(1)]); % Periodic.
   
   H(n) = dx * sum(1/2 * p.^2 + 1/2 * qx.^2 + 1 - cos(q));
end