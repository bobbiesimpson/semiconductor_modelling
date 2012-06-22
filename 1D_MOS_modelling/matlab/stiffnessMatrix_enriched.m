function [ K_el ] = stiffnessMatrix_enriched(x1,x2,eps_1,shifting,xi,eps_2)
%   Evaluation of enriched stiffness matrix
%   The input variables are as follows:
%   x1, x2      - Element coordinates
%   eps_1       - Element permittivity
%   shifting    - BOOL value to determine whether nodal shifting is used
%   xi          - X-coord of interface
%   eps_2       - (Optional) Is needed when xi lies in element

K_el = zeros(2,2);
l = x2-x1;               % length of element
ngp = 2;
[gpt gwt] = gaussQuad(ngp);
dN1 = -1/l;
dN2 = 1/l;
    
if(xi < x2 && xi > x1)        %   Interface lies in mid-element

    l_sub_mat = [xi-x1 x2-xi];      %   Sub element lengths
    b = [xi x2];
    a = [x1 xi];
    eps = [eps_1 eps_2];
    
    for sub_el=1:2
        lsub = l_sub_mat(sub_el);
        for i=1:ngp
            %   coord transformation
            [N1 N2] = linearShapeFn(gpt(i));
            x = [N1 N2] * [a(sub_el) b(sub_el)]';
            N1 = -(x - x2)/l;
            N2 = (x - x1)/l;
            if(shifting) 
                enrichTerm1 = hatEnrichment(x,xi,x1);
                enrichTerm2 = hatEnrichment(x,xi,x2);
            else
                enrichTerm1 = hatEnrichment(x,xi);
                enrichTerm2 = enrichTerm1;
            end
            deriv1 =  N1*hatEnrichmentDeriv(x,xi) + ...
                             dN1*enrichTerm1;
            deriv2 =  N2*hatEnrichmentDeriv(x,xi) + ...
                             dN2*enrichTerm2;
                         
            K_el(1,1) = K_el(1,1) + deriv1*eps(sub_el)*deriv1*lsub/2*gwt(i); 
            K_el(1,2) = K_el(1,2) + deriv1*eps(sub_el)*deriv2*lsub/2*gwt(i);
            K_el(2,1) = K_el(1,2);
            K_el(2,2) = K_el(2,2) + deriv2*eps(sub_el)*deriv2*lsub/2*gwt(i);
        end
    end

else
    %   no sub-elements
    for i=1:ngp
        [N1 N2] = linearShapeFn(gpt(i));
        x = [N1 N2] * [x1 x2]';
        if(shifting)
            enrichTerm1 = hatEnrichment(x,xi,x1);
            enrichTerm2 = hatEnrichment(x,xi,x2);
        else
            enrichTerm1 = hatEnrichment(x,xi);
            enrichTerm2 = enrichTerm1;
        end
        deriv1 =  N1*hatEnrichmentDeriv(x,xi) + ...
            dN1*enrichTerm1;
        deriv2 =  N2*hatEnrichmentDeriv(x,xi) + ...
            dN2*enrichTerm2;
        
        K_el(1,1) = K_el(1,1) + deriv1*eps_1*deriv1*l/2*gwt(i);
        K_el(1,2) = K_el(1,2) + deriv1*eps_1*deriv2*l/2*gwt(i);
        K_el(2,1) = K_el(1,2);
        K_el(2,2) = K_el(2,2) + deriv2*eps_1*deriv2*l/2*gwt(i);
    end
end


