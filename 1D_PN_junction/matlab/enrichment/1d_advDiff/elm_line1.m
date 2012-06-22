function xe = elm_line1(x1, x2, n, ratio)

% here we are doing some mesh grading where we have the following givens
%
%   - Length (between the start and finish points), given by x2 - x1
%   - number of elements we want (n)
%   - ratio between last element length and first element length

%  what we return is xe: element end nodes

if(n==1)    % we have one element
    xe(1) = x1;
    xe(2) = x2;
    return;
end

% *************************************************
% ***** we get here if we have > 1 element ********
%**************************************************

if(ratio==1)
    alpha = 1.0;
    factor = 1.0 / n;
else
    exp = 1 / (n-1);
    alpha = ratio^exp;
    factor = ( 1 - alpha ) / ( 1- alpha^n);
end

deltax = (x2 - x1) * factor;
xe(1) = x1;

for(k=2:n+1)
    xe(k) = xe(k-1) + deltax;
    deltax = deltax * alpha;
end

return;