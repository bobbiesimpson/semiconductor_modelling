function [ chi ] = hatEnrichment( x, x_i, nodalX )
% This is our hat enrichment function with the following arguments
%   x           Our 1D x-coordinate
%   xi          The c-coordinate of our interface
%   nodalX      If we are using nodal shifting, the x-coord of the relevant
%               node

%   We can use a variable number of arguments in this function where if we
%   leave out the nodalX value we do not use nodal shifting

if (nargin == 2)
    chi = abs( x - x_i);
elseif (nargin == 3)
    nodalChiValue = abs(nodalX - x_i);
    chi = abs( x - x_i) - nodalChiValue;
end

