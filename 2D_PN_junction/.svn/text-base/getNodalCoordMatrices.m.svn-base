function [ xcoordM, ycoordM, nodalNumbering ] = getNodalCoordMatrices( nodes, elConn, nex_p, nex_n,  ney_t )
    %   pass in the nodal coordindates and element connectivity and form
    %   the matrices of the x and y coordinates which are required for the
    %   'surf' function

    % p region
    n_els_p=nex_p*ney_t;    % num elements in p-region
    n_els_n=nex_n*ney_t;    % num elements in n-region
    
    xcoordM=zeros(ney_t+1,nex_p+nex_n+1);   % where we store the x-coordinates
    ycoordM=zeros(ney_t+1,nex_p+nex_n+1);   % and likewise for the y-coords
    nodalNumbering=zeros((ney_t+1)*(nex_n+nex_p+1),1);  % this will contain a vector of the nodes in the order they appear in the coord matrices
    
    colCounter=1;
    
    for col=1:nex_p
        elements=col:nex_p:n_els_p;
        colNodalNums=[elConn(elements,4)' elConn(elements(size(elements,2)),3)'];
        
        %colNodalNums=unique(elConn(col:nex_p:n_els_p,3:4))
        row=colCounter*(ney_t+1) - ( ney_t + 1 ) + 1;
        nodalNumbering(row:row+ney_t)=colNodalNums;
        xcoordM(:,col)=nodes(colNodalNums,1);
        ycoordM(:,col)=nodes(colNodalNums,2);
        
        colCounter=colCounter+1;
    end
    
    elements=n_els_p+1:nex_n:n_els_n+n_els_p;
    colNodalNums=[elConn(elements,4)' elConn(elements(size(elements,2)),3)']; 
    
    %colNodalNums=unique(elConn(n_els_p+1:nex_n:n_els_n+n_els_p,3:4));
    row=colCounter*(ney_t+1) - ( ney_t + 1 ) + 1;
    nodalNumbering(row:row+ney_t)=colNodalNums;
    xcoordM(:,nex_p+1)=nodes(colNodalNums,1);
    ycoordM(:,nex_p+1)=nodes(colNodalNums,2);
    
    colCounter=colCounter+1;
        
    for col=1:nex_n
        elements=col+n_els_p:nex_n:n_els_n+n_els_p;
        colNodalNums=[elConn(elements,4)' elConn(elements(size(elements,2)),3)']; 
        
        %colNodalNums=unique(elConn(col+n_els_p:nex_n:n_els_n+n_els_p,1:2));
         row=colCounter*(ney_t+1) - ( ney_t + 1 ) + 1;
        nodalNumbering(row:row+ney_t)=colNodalNums;
        xcoordM(:,col+nex_p+1)=nodes(colNodalNums,1);
        ycoordM(:,col+nex_p+1)=nodes(colNodalNums,2);
        
        colCounter=colCounter+1;
    end

end

