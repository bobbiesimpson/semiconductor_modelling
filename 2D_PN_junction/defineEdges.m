function [ leftedge, bottomedge, topedge, rightedge, junction_edge, corners ] = defineEdges(  nex_t, ney_t, nex, ney, nodeNum_p, nodeNum_n, nodes)
%Define the edges of the 2D pn junction mesh

uln=1;                      %   upper left node
urn=(nex_t)*(ney_t+1)+1;    %   upper right node
lln=ney_t+1;                %   lower left node
lrn=length(nodes);          %   lower right node
corners=[lln lrn urn uln];

leftedge=1:ney_t+1;         %   nodes on left edge etc
bottomedge=ney_t+1:ney_t+1:length(nodes);
topedge=1:ney_t+1:(nex_t+1)*(ney_t+1);
rightedge=(ney_t+1)*(nex_t)+1:length(nodes);
junction_edge=(ney+1)*(nex)+1:nodeNum_p;

end

