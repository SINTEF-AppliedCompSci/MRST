function faces=findfracboundaryfaces(G,tol)
% Generates a list of global faces belonging to fractures. These global
% faces intersect with the external boundary with the domain. This is to
% facilitate assigning boundary conditions to these faces. The list of
% faces can be plugged into addBC. The inputs are the global grid which
% needs to contain fracture grid data.
%
% The output is a struct class object containing the fields West, East,
% South, North, Top and Bottom, with each field containing a list of
% fracture faces.
%
% Note that this function only works for 3D rectangular domains.

faces=struct('West',[],'East',[],'South',[],'North',[],'Top',[],'Bottom',[]);

%% Extract external faces
mincoords=min(G.Matrix.nodes.coords);
maxcoords=max(G.Matrix.nodes.coords);

% West face
westpoint=mincoords;
westnormal=[1 0 0];

% East face
eastpoint=maxcoords;
eastnormal=[-1 0 0];

% South face
southpoint=mincoords;
southnormal=[0 1 0];

% North face
northpoint=maxcoords;
northnormal=[0 -1 0];

% Top face
toppoint=mincoords;
topnormal=[0 0 1];

% Bottom face
botpoint=maxcoords;
botnormal=[0 0 -1];


%% Extract boundary faces
numfracs=length(fieldnames(G.FracGrid));
for i=1:numfracs
    fieldname=['Frac',num2str(i)];
    F=G.FracGrid.(fieldname);
        
    % Check if the nodes lie on the external boundary of the domain. If so,
    % mark which boundary they are on. Boundary tags are from columns 1:6,
    % sequentially representing west, east, south, north, top, bottom. Note
    % that during fracture grid generation, the flat fracture plane was
    % extruded in both directions of the fracture by aperture/2. This
    % process is reversed before checking if the nodes lie on the external
    % boundary. 
    numnodes=F.nodes.num/2;
    nodes=0.5*(F.nodes.coords(1:numnodes,:)+...
        F.nodes.coords(numnodes+(1:numnodes),:));
    xsectflag=false(numnodes,6);
    for j=1:numnodes
        nodej=nodes(j,:);
        
        % West face check
        xsectflag(j,1)=pointplanedistance(nodej,westnormal,westpoint)<tol;
        
        % East face check
        xsectflag(j,2)=pointplanedistance(nodej,eastnormal,eastpoint)<tol;
        
        % South face check
        xsectflag(j,3)=pointplanedistance(nodej,southnormal,southpoint)<tol;
        
        % North face check
        xsectflag(j,4)=pointplanedistance(nodej,northnormal,northpoint)<tol;
        
        % Top face check
        xsectflag(j,5)=pointplanedistance(nodej,topnormal,toppoint)<tol;
        
        % Bottom face check
        xsectflag(j,6)=pointplanedistance(nodej,botnormal,botpoint)<tol;
    end
    xsectflag=repmat(xsectflag,2,1); % 'expand' back
    
    
    numfaces=F.faces.num;
    facestart=F.faces.start;
    for j=1:numfaces
        % Skip if face is internal in the fracture (this is satisfied if
        % all face neighbours are nonzero)
        if all(F.faces.neighbors(j,:))
            continue;
        end
        
        % Extract node indices of a face
        npos=F.faces.nodePos(j):(F.faces.nodePos(j+1)-1);
        nodeind=F.faces.nodes(npos);
        
        % If all nodes intersect a boundary, face intersects the boundary
        if all(xsectflag(nodeind,1)) % West
            faces.West(end+1)=facestart-1+j;
        elseif all(xsectflag(nodeind,2)) % East
            faces.East(end+1)=facestart-1+j;
        elseif all(xsectflag(nodeind,3)) % South
            faces.South(end+1)=facestart-1+j;
        elseif all(xsectflag(nodeind,4)) % North
            faces.North(end+1)=facestart-1+j;
        elseif all(xsectflag(nodeind,5)) % Top
            faces.Top(end+1)=facestart-1+j;
        elseif all(xsectflag(nodeind,6)) % Bottom
            faces.Bottom(end+1)=facestart-1+j;
        end
    end
   
end


end