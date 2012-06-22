function xe = twoGradedMesh(x1,x2,x3,neL,lRatio,neR, rRatio);

%===============================================
%   This meshes a 1D junction with the following parameters
%   x1, x2, x3:   Coordinates of left contact, junction, and right contact
%   ne_l          No. elements in left region
%   lRatio        Ratio of last element/first element in left region
%   ne_r          No. elements in right region
%   rRatio        Ratio of first element/last element in right region
%===============================================

%------------
% left region
%------------
xe(1) = x1;
lRatio = 1/lRatio;  %   invert ratio for LHS (grade towards junction)
if(neL==1)
  xe(2) = x2;
else
    if(lRatio==1)
        alpha = 1.0;
        factor = 1.0/neL;
    else
        texp = 1/(neL - 1);
        alpha = lRatio^texp;
        factor = (1.0-alpha)/(1.0-alpha^neL);
    end
    
    deltax=(x2-xe(1))*factor;   % length of the first element
    for i=2:neL+1
        xe(i)  = xe(i-1)+deltax;
        deltax = deltax*alpha;
    end
end

%--------------
% right region
%--------------
stNod = neL+1;
if(neR==1)
    xe(stNod+1) = x3;
    return;
end

if(rRatio==1)
  alpha  = 1.0;
  factor = 1.0/neR;
else
  texp = 1/(neR-1);
  alpha = rRatio^texp;
  factor = (1.0-alpha)/(1.0-alpha^neR);
end

deltax=(x3-xe(stNod))*factor;   % length of the first element

for i=stNod+1:neR+stNod
  xe(i)  = xe(i-1)+deltax;
  deltax = deltax*alpha;
end

%-----
% done
%-----

return;