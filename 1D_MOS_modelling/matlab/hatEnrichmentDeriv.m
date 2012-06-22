function [ chiDeriv ] = hatEnrichmentDeriv( x, x_i )
% This is our hat enrichment function derivative with the following arguments
%   x           Our 1D x-coordinate
%   xi          The c-coordinate of our interface

%   We can use a heaviside function to define the gradient of the hat
%   function. BUT be warned that errors will occur if gauss points lie
%   exactly on the interface
chiDeriv = 2*heaviside(x - x_i)-1;
