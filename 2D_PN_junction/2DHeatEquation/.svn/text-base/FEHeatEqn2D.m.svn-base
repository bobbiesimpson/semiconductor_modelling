% To test my coding abilites in 2D this scipt will solve the heat equation
% on a square domain with Dirichlet boundary conditions on the lower left
% sides and Neumann conditions on the upper right.

%           q.n=2
%       -------------
%       |           |
%       |           |
%  T=10 |           | q.n=2
%       |           |
%       -------------
%           T=10

% -------------------------------
% ----- Problem setup -----------
% -------------------------------
lx=10; ly=10;
k=5;                                 % conductivity
D=[k 0; 0 k];                        % D matrix of conductivity (isotropic)
Timposed=10;
fluxImposed=2;

% -------------------------------
% ----- Generate the mesh -------
% -------------------------------
nx=100; ny=100;                         %   number of elements in each direction
[NodeNum, ElemNum, xx, yy, elConn] = GetMesh(lx, ly, nx, ny);
coords=[xx yy]

% ----------------------------------------
% ----- Define the Neumann boundaries ----
% ----------------------------------------
uln=1;                              % upper left node            
lln=ny+1;                           % lower left node
urn=(ny+1)*(nx+1)-ny;               % upper right node
lrn=length(nodes);                  % lower right node

topedge=[uln:ny+1:urn];             % nodes on topedge
rightedge=[urn:lrn];                % nodes on right edge

% -------------------------------------
% ----- Define the dirichlet nodes ----
% -------------------------------------
leftedge=[1:ny+1];                         % left edge nodes
bottomedge=[ny+1:ny+1:(ny+1)*(nx+1)];      % bottom edge nodes
fixedNodes=unique([leftedge bottomedge]);   % the nodes on dirichlet boundaries
Tfixed=zeros(size(fixedNodes));            % set the Temp. to equal zero

% ----------
% plot the mesh
% ----------
PlotMesh(elConn,xx,yy,ElemNum)

% ---------------------
% ----- processing ----
% ---------------------
d=zeros(NodeNum,1);
f=zeros(NodeNum,1);
K=sparse(NodeNum,NodeNum);

% compute boundary force vector



for el=1:ElemNum
    
end








