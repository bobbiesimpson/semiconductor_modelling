function [ L2norm ] = L2norm( nE, meshcoords, psi, enr_el_conn, el_conn, shifting)
%   Function to calculate the L2 norm of our current Psi solution
%   Works with both enriched and unenriched solutions

load('mesh_ref.mat');
load('psi_ref.mat');

meshref = mesh_ref;
psiref = psi_ref;

global q eps Na Nd x_jun nep Vbi xn xp

L2norm = 0;  

for element = 1:nE
    loc_nod_enr = find(enr_el_conn(element,:));          % local enriched nodes #s
    glb_dof_enr = enr_el_conn(element,loc_nod_enr);      % global enriched nod #s
    
    nodes = el_conn(element,:);
    lcoords = meshcoords(nodes);                         % local node coordinates
    el_length = lcoords(2) - lcoords(1);
    detJacob = el_length / 2;
    
    if(element < nep + 1)
        doping = -Na;                                    %   p-type doping
    else
        doping = Nd;                                     %   n-type doping
    end
    
    %   get the nodal values of the enrichment functions 
    if ~isempty(loc_nod_enr) && shifting
        enrichFn_nod = zeros(1,length(nodes));
        enrichFnDeriv_nod = zeros(1,length(nodes));
        for n = 1:length(nodes)
            [enrichFn_node(n) enrichFnDeriv_node(n)] = pn_enrichment_function(lcoords(n), x_jun, xp, xn, eps, doping, Vbi, q);
        end
    end
    
    if isempty(loc_nod_enr) 
      [gaussPoint gaussWeight] = gaussQuad(2);   %   get gauss pts and weights                                      %   gauss pts for unenriched elmt
    else
      [gaussPoint gaussWeight] = gaussQuad(10);  %   get gauss pts and weights                                     %   gauss pts for enriched elmt
    end

    for gp = 1:length(gaussPoint)
        gpt = gaussPoint(gp);
        gwt = gaussWeight(gp);
        shape = linearShapeFn(gpt);
        xcoord = shape * lcoords';                  %   coord of GP
        interp_psi = shape * psi(nodes(:));         %   unenriched interpolated psi
        
        enrichFn = 0;
        enrichFnDeriv = 0;
        if ~isempty(loc_nod_enr)                    % calculate the enrichment functions and derivatives
            [enrichFn enrichFnDeriv] = pn_enrichment_function(xcoord, x_jun, xp, xn, eps, doping, Vbi, q);
            
            if shifting
                enrichFn = [enrichFn - enrichFn_node(1) enrichFn - enrichFn_node(2)];
                enrichFnDeriv = [enrichFnDeriv - enrichFnDeriv_node(1) enrichFnDeriv - enrichFnDeriv_node(2)];
                interp_psi = interp_psi + shape(loc_nod_enr) .* enrichFn(loc_nod_enr) ...
                    * psi(glb_dof_enr);
            else
                interp_psi = interp_psi + shape(loc_nod_enr).*enrichFn * currentPsi(glb_dof_enr);
            end
        end
        
        %   we need to interpolate the reference solution - first job is to
        %   find the nodes of the element we are on

        [val index1] = min(abs(meshref - xcoord));
        temp = meshref(index1);
        meshref(index1) = 1e20;                         %   put arbitrary large number in place
        [val index2] = min(abs(meshref - xcoord));
        meshref(index1) = temp;
        
        loc_nod_ref = sort([index1 index2]);            %   the global node #s of our element
        lcoords_ref = meshref(loc_nod_ref);            %   the coordinates of our element
        
        xi = 2 * (xcoord - lcoords_ref(1)) / (lcoords_ref(2) - lcoords_ref(1)) - 1;
        shapeRef = linearShapeFn(xi);
        psi_nod_ref = psiref(loc_nod_ref);
        psi_ref_gp = shapeRef * psi_nod_ref;
        
        L2norm = L2norm + (interp_psi-psi_ref_gp)^2 * gwt * detJacob;
    end
    
end

L2norm = sqrt(L2norm);

end

