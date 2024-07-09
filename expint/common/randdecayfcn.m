function varargout = randdecayfcn(regularity, modes)
% RANDDECAYFCN - Returns a random vector with  
%             prescribed regularity, or more precise; 
%             with a prescribed decay of the Fourier coefficients.
%             *Vector returned is in Fourier space*.
%
% SYNOPSIS:
%   vector           = randdecayfcn(regularity, modes)
%   [vector, string] = randdecayfcn(regulatity, modes)
%
% PARAMETERS:
%   regularity  - Desired regularity (polynomial order) of returned
%                 vector. This could be either a string in the format
%                 'reg3', or just the integer 3.
%   modes       - Number of Fourier modes to include in returned vector.
%
% RETURNS:
%   vector - A vector in Fourier space, and with the given
%           decay rate. Normalised such that the norm will be 1 after
%           an IFFT.
%   regstring - A character-representation of the regularity number.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.3 $ $Date: 2005/10/12 16:28:00 $

if isa(regularity, 'char'),
   [start, finish, tok] = regexp(lower(regularity), 'reg\D*0*(\d+)');
   regstring = regularity(tok{1}(1) : tok{1}(2));         % Highly convoluted...
   regularity = str2num(regstring);
elseif isa(regularity, 'numeric'),
   % no-op.
else
   error('randdecayfcn:badopt', 'Unknown type of regularity argument');
end

m  = modes / 2;

% Complex vector with frequencies in fourier domain.
Di = -1i * [0, 1./(1:m-1), 0, 1./(-m+1:-1)]';

% Generate random vector and enforce correct decay rate:
v = (Di.^regularity) .* (randn([modes, 1]) + i*randn([modes, 1]));

% The vector is normalized such that it will be normalized
% after an ifft, therefore the sqrt-factor.
v = sqrt(modes) * v / norm(v);

varargout{1} = v;

if nargout == 2,
   varargout{2} = num2str(regularity);
end
