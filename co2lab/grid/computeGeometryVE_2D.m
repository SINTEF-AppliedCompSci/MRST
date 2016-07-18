function Gt = computeGeometryVE_2D(Gt, varargin)
%Compute geometry of grid of top-surface grid so that the geometry is
% defined by its projection and the z value is the value at face
% corresponting to the top of the column
% 
%
% SYNOPSIS:
%   Gt = computeGeometryVE(Gt)
%   Gt = computeGeometryVE(Gt, varargin)
%
% PARAMETERS:
%   Gt       - Grid structure as defined by function 'topSurfaceGrid'.
%              Gt.parent is the original grid
%   
%
%   'pn'/pv - List of 'key'/value pairs for supplying optional parameters.
%             The supported options are
%               - Verbose -- Whether or not to display verbose output as
%                       the process progresses.  Possible values are
%                       TRUE and FALSE.  Default value equals mrstVerbose.
%
% RETURNS:
%   Gt - Grid structure with added fields:
%         - cells
%             - volumes
%             - centroids
%             - z
%         - faces
%             - areas
%             - normals                !!! AREA WEIGHTED !!!
%             - centroids
%             - z

%
% COMMENTS:
%   PLEASE NOTE: Face normals have length equal to face areas

%{
#COPYRIGHT#
%}

% $Date: 2012-01-30 11:41:03 +0100 (Mon, 30 Jan 2012) $
% $Revision: 9020 $

% Setup
assert(size(Gt.faces.nodes, 2)==1);
assert(Gt.griddim==2)
assert(numel(Gt.columns.cells)==Gt.parent.cells.num)

% compute2D geometry
type = Gt.type;
Gt=computeGeometry(Gt);
Gt.cells.z= Gt.parent.faces.centroids(Gt.cells.map3DFace, 3);
%faceCentroids = (coords(faceEdges(:,2),:)+ coords(faceEdges(:,1),:))/2;
%Gt.faces.z=faceCentroids(:,3);

faceEdges = reshape(Gt.faces.nodes,2,[])';
z = (Gt.nodes.z(faceEdges(:,2))+ Gt.nodes.z(faceEdges(:,1)))/2;
Gt.faces.z=z;
%Gt.type = [Gt.parent.type, mfilename];
Gt.type = [type, mfilename]; % @@ Was there a reason for refering to the
                                % parent type instead?  
