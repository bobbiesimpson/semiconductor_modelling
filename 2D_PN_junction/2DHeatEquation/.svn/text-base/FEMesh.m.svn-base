function [node,element,typE] = FEMesh()

% filename: FEMesh.m
% purpose: read FE mesh data. meshing done with GID.
%
% Author: S Natarajan
% November 2010
% -------------------------------------------------------------------------

% open file with nodal coordinate information. get the number of lines,
% this will be the total number of nodes in the system...
% first column -> node number, column 2 -4 -> nodal coordinate X,Y,Z. we
% have set ndim = 2, so it only reads X,Y coordinates.
fid = fopen('nodes.inp','r') ;
count = 0 ;
while ~feof(fid)
    line = fgetl(fid) ;
    if isempty(line) 
        continue
    end
    count = count + 1 ;
end

numnode = count ;
ndim = 2 ;
fclose(fid) ;

fid = fopen('nodes.inp','r') ;
for i=1:numnode
    tmp = str2num(fgets(fid)) ;
    [n, node(n,:)] = deal(tmp(1),tmp(2:1+ndim)) ;
end

fclose(fid) ;

% read element data. open file with element connectivity and corresponding
% material id associated with the element. this will be the last coulmn.
% first column -> element number, column 2 - 5 -> element connectivity,
% column 6 -> material id.
fid = fopen('elms.inp','r') ;

count = 0 ;
while ~feof(fid)
    line = fgetl(fid) ;
    if isempty(line) 
        continue
    end
    count = count + 1 ;
end

numelem = count ;
numcol = 4 ;
fclose(fid) ;

fid = fopen('elms.inp','r') ;

kold = 0;
kold = 0;
for i=1:numelem
tmp = str2num(fgets(fid));
knew = kold + numcol;
[ncell,element(ncell,:),typE(ncell,1)] = deal(tmp(1),...
    tmp(2:1+numcol),tmp(2+numcol));
kold = knew;
end

fclose(fid) ;