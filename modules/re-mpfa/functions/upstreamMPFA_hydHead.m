function [krwUp] = upstreamMPFA_hydHead(krw,F,boundF,G_struct,bc_struct,bc_val,hydHead,zetaCells,zetaFaces)

% Computes the upstream weighting of the relative permeability 
%
% SYNOPSIS:
%   krwAr = upstreamMPFA_hydHead(krw,gradF,boundF,G_struct,bc_struct,bc_val,hydHead,zetaCells,zetaFaces)
%
% PARAMETERS:
%   krw         - Function, relative permeability function krw = krw(psi)
%   F           - Function, MPFA discrete operator that acts on the
%                 hydraulic head
%   boundF      - Function, MPFA discrete operator that acts on bc_val 
%   G_struct    - Structure, Grid structure from MRST
%   bc_struct   - Structure, Boundary conditions structure from MRST
%   bc_val      - Vector, containing the values of boundary conditions.
%                 This vector must have a length equal to G_struct.faces.num            
%   hydHead     - Vector, containing the values of hydraulic head. This
%                 vector must have a length equal to G_struct.cells.num
%   zetaCells   - Vector, containing the values of the elevation head at
%                 the cell centers. This vector must have a length equal to
%                 G_struct.cells.num
%   zetaFaces   - Vector, contaning the values of the elevation head at
%                 the faces centroids. This vector must have a length equal
%                 to G_struct.faces.num
%
%  RETURNS:
%   krwUp       - Vector, containing the upstream weighted relative
%                 permeabilities at the faces. This vector will have a
%                 lenght equal to G_struct.faces.num
%


%% Exctracting grid and boundary information
fNei = G_struct.faces.neighbors;               % extracting faces neighbors
int_fNei = fNei(all(fNei ~= 0,2),:);    % internal faces neighbors
int_f = find(all(fNei ~= 0,2));         % internal faces
neumFaces = bc_struct.face(all(bc_struct.face.*strcmp(bc_struct.type','flux') ~= 0,2),:); % Extracting the indices of the Neumann faces
diriFaces = bc_struct.face(all(bc_struct.face.*strcmp(bc_struct.type','pressure') ~= 0,2),:); % Extracting the indices of the Dirichlet faces
krwUp = zeros(G_struct.faces.num,1);           % initializing the upstream weighted krw

%% Neumann boundaries relative permeabilities
krwUp(neumFaces) = 1;   % Since the flux is already known, we set krw=1

%% Computing the flux at every face of the grid
flux = F(hydHead) + boundF(bc_val);

%% Dirichlet boundaries relative permeabilities
hFacesIndex = find(ismember(bc_struct.face,diriFaces));  % extracting h at every dirichlet face from the bc struct
diriNeigh = fNei(diriFaces,:)'; % neighboring cells of dirichlet faces
diriCells = diriNeigh(diriNeigh > 0); % cells corresponding to each dirichlet face
upIsFace = fNei(diriFaces,1) == 0;   % is the upstream direction a face? ~upIsFace corresponds to a cell
downIsFace = fNei(diriFaces,2) == 0; % is the downstream direction a face? ~downIsFace corresponds to a cell
% This may look complicated, but is the easiest way to do it. Basically we
% compute the flux at each dirichlet face. If the flux >= 0 then we must
% evaluate in the upstream direction. Then, we ask if the upstream
% direction is a face or not. If true, then we evaluate at the face, if not
% we evaluate at the cell center. If the flux < 0, then we must evaluate
% in the downstream direction. Then, we ask if the downstream direction is 
% a face or not. If true, then we evaluate at the face, if not we evaluate
% at the cell center
krwUp(diriFaces) = ( ...
    (flux(diriFaces) >= 0) .* upIsFace .* krw(bc_struct.value(hFacesIndex) - zetaFaces(diriFaces)) + ...
    (flux(diriFaces) >= 0) .* ~upIsFace .* krw(hydHead(diriCells) - zetaCells(diriCells)) + ...
    (flux(diriFaces) < 0) .* downIsFace .*  krw(bc_struct.value(hFacesIndex) - zetaFaces(diriFaces)) + ...
    (flux(diriFaces) < 0) .* ~downIsFace .* krw(hydHead(diriCells) - zetaCells(diriCells)) ...
    );
%% Internal faces
% upCenters determines the indices at which the relative permeabilities must
% be evaluated. If flux >= 0 then krw is evaluated at the upstream
% direction. If flux < 0 then krw is evaluated at the downstream direction.
upCenters = ( ...
    (flux(int_f) >= 0) .* int_fNei(:,1) + ... % evaluate at i-1
    (flux(int_f) < 0)  .* int_fNei(:,2)   ... % evaluate at i
    );
krwUp(int_f) = krw(hydHead(upCenters) - zetaCells(upCenters)); % here, we compute the upstream krw for every internal face
end
