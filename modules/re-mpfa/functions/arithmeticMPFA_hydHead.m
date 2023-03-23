function [krwAr] = arithmeticMPFA_hydHead(krw,G_struct,bc_struct,hydHead,zetaCells,zetaFaces)
% Computes the arithmetic average of the relative permeability 
%
% SYNOPSIS:
%   krwAr = arithmeticMPFA_hydHead(krw,G_struct,bc_struct,hydHead,zetaCells,zetaFaces)
%
% PARAMETERS:
%   krw         - Function, relative permeability function krw = krw(psi)
%   G_struct    - Structure, Grid structure from MRST
%   bc_struct   - Structure, Boundary conditions structure from MRST
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
%   krwAr       - Vector, containing the arithmetic averaged relative
%                 permeabilities at the faces. This vector will have a
%                 lenght equal to G_struct.faces.num
%



fNei = G_struct.faces.neighbors;        % extracting faces neighbors
int_fNei = fNei(all(fNei ~= 0,2),:);    % internal faces neighbors
int_f = find(all(fNei ~= 0,2));         % internal faces
neumFaces = bc_struct.face(all(bc_struct.face.*strcmp(bc_struct.type','flux') ~= 0,2),:); % extracting neumann faces
diriFaces = bc_struct.face(all(bc_struct.face.*strcmp(bc_struct.type','pressure') ~= 0,2),:); % extracting dirichlet faces

krwAr = zeros(G_struct.faces.num,1);

% Neumann boundaries relative permeabilities
krwAr(neumFaces) = 1; % note that we don't really use this values in the computation since the fluxes are known

% Dirichlet boundaries relative permeabilities
hFacesIndex = find(ismember(bc_struct.face,diriFaces));  % extracting index of every dirichlet face from the bc struct
diriNeigh = fNei(diriFaces,:)'; % neighboring cells of dirichlet faces
diriCells = diriNeigh(diriNeigh > 0); % cells corresponding to each dirichlet face
krwAr(diriFaces) = 0.5 .* ( ...
                                krw(bc_struct.value(hFacesIndex) - zetaFaces(diriFaces)) + ...
                                krw(hydHead(diriCells) - zetaCells(diriCells)) ...
                           );
% Internal faces
krwAr(int_f) = 0.5 .* (...
                          krw(hydHead(int_fNei(:,1))-zetaCells(int_fNei(:,1))) + ...
                          krw(hydHead(int_fNei(:,2))-zetaCells(int_fNei(:,2))) ...
                       );                                                  
end
