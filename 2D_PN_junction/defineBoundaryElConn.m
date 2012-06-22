function [ leftEdgeConn, bottomEdgeConn, rightEdgeConn, topEdgeConn ] = defineBoundaryElConn(corners, nex_t, ney_t )
%On the boundary, we need the element connectivity matrices for each of the
%edge to allow the flux vector to be integrated

lln=corners(1); lrn=corners(2); urn=corners(3); uln=corners(4);
rightEdgeConn=[urn:lrn-1; urn+1:lrn]';      % connectivity matrix of right edge
topEdgeConn=[uln:ney_t+1:urn-(ney_t+1); uln+ney_t+1:(ney_t+1):urn]';    % same for top
bottomEdgeConn=[lln:ney_t+1:lrn-(ney_t+1); lln+ney_t+1:ney_t+1:lrn ]';
leftEdgeConn=[1:lln-1; 2:lln]';


end

