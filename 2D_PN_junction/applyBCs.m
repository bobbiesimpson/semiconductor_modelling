function [ K, f ] = applyBCs( K, f, fixedNodes, fixedNodeValues )
% Given a matrix K and the associated body foce vector f, we apply the
% boundary conditions and return the modified matrix and vector

bcwt = trace(K)/length(K);
K(fixedNodes,:) = 0;
f = f- K(:,fixedNodes) * fixedNodeValues;
K(:,fixedNodes) = 0;
K(fixedNodes,fixedNodes) = speye(length(fixedNodes))*bcwt;
f(fixedNodes) = fixedNodeValues*bcwt;

end

