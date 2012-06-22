function xe = gradedMesh(x1,x2,oxthick,nox,n,ratio);

%===============================================
% Unsymmetrical discretization of a line segment
% into a graded mesh of n elements
% subtended between the left end-point x1
% and the right end-point x2
% oxthick is the thickness of the oxide
%
% The variable "ratio" is the ratio
% of the length of the last
% element to the length of the first element.
%
% The number of elements on the oxide layer is nox
%
% xe: element end-nodes
%===============================================

%------------
% oxide layer -> assume evenly spaced over oxide
%------------

if(nox==1)
  xe(1) = x1;
  xe(2) = x1 + oxthick;
else
  xe(1) = x1;
  step = oxthick / nox;
  for i=2:(nox+1) 
    xe(i) = xe(i-1) + step;
  end
end

%--------------
% semiconductor
%--------------

stNod = nox + 1;
if(n==1)
    xe(stNod+1) = x2;
    return;
end

if(ratio==1)
  alpha  = 1.0;
  factor = 1.0/n;
else
  texp = 1/(n-1);
  alpha = ratio^texp;
  factor = (1.0-alpha)/(1.0-alpha^n);
end

deltax=(x2-xe(stNod))*factor;   % length of the first element

for i=stNod+1:n+stNod
  xe(i)  = xe(i-1)+deltax;
  deltax = deltax*alpha;
end

%-----
% done
%-----

return;