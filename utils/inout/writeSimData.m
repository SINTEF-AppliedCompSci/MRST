function writeSimData(varargin)
%Write SAM simulator grid input files.
%
% SYNOPSIS:
%   writeSimData(G, caseName)
%   writeSimData(G, caseName, 'src', src, 'rock', rock, ...)
%
% PARAMETERS:
%   G - Grid structure as detailed in grid_structure
%   caseName - Simulation case base name.  String assumed to contain a
%              (relative) base name.
%
%              files:
%                 Geometry: [caseName, '-geom.dat']
%                 Topology: [caseName, '-topo.dat']
%                 %IJK-map:  [caseName, '-map.dat']
%
% RETURNS:
%   Nothing
%
% SEE ALSO:
%   `fopen`, `fprintf`.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

   opt = struct('casename', 'case', ...
                'grid', [], ...
                'rock', [], ...
                'src', [],...
                'wells',[]);
   opt = merge_options(opt, varargin{:});

   pstr = fileparts(opt.casename);
   if ~isempty(pstr) && ~(exist(pstr, 'dir')==7),
      mkdir(pstr);
   end
   if ~isempty(opt.grid), writeSimGrid(opt.casename, opt.grid); end
   if ~isempty(opt.rock), writeSimRock(opt.casename, opt.rock); end
   %if ~isempty(opt.wells), writeSimWells(opt.casename, opt.wells); end
   if (~isempty(opt.src) || ~isempty(opt.wells)), writeSimSrcAndWells(opt.casename, opt.src,opt.wells); end

end
function writeSimSrcAndWells(caseName, src, W)
   fp = fopen(sprintf('%s-wells.dat', caseName),'w');
   if(~fp)
      error('Could not open file')
   end
   num_wells=0;
   if(~isempty(src))
      src.rate = convertTo(src.rate, meter^3/day);
      num_wells = num_wells + numel(src.cell);
   end
   if(~isempty(W))
      num_wells = num_wells + numel(W);
   end
   fprintf(fp, '%d\n', num_wells);
   if(~isempty(src))
      for i=1:numel(src.cell),
         if src.rate(i)>0,
            fprintf(fp, 'SimpleInjection I%d 1 %15.15g %15.15g; %5d %15.15g\n', ...
               i, src.rate(i), src.sat(i), src.cell(i)-1, src.rate(i));
         else
            fprintf(fp, 'SimpleProduction P%d 1 %15.15g %15.15g; %5d %15.15g\n', ...
               i, -src.rate(i), src.sat(i), src.cell(i)-1, -src.rate(i));
         end
      end
   end
   if(~isempty(W))
      for i=1:numel(W),
         w=W(i);
         well_type=w.type;
         switch lower(well_type),
            case {'bhp'}
               fprintf(fp, 'Pressure %s %d %15.15g %4.4g',w.name,...
                  numel(w.cells),w.val,w.compi(1));
            case {'rate'}
               if(w.val>0)
                  fprintf(fp, 'Injection %s %d %15.15g %4.4g',w.name,...
                     numel(w.cells),convertTo(w.val,meter^3/day),w.compi(1));
               else
                  fprintf(fp, 'Production %s %d %15.15g %4.4g',w.name,...
                     numel(w.cells),abs(convertTo(w.val,meter^3/day)),w.compi(1));
               end
            otherwise
               error('Wrong well');
         end
         wunit=centi*poise * meter^3 / (day * barsa);
         for i=1:numel(w.cells),
            fprintf(fp,'; %15.15g %15.15g',w.cells(i)-1,convertTo(w.WI(i),wunit));
         end
         fprintf(fp,'\n');
      end
   end
end
function writeSimRock(caseName, rock)
fp = fopen(sprintf('%s-perm.dat', caseName),'w');
if(fp)
   fprintf(fp, '%d\n', size(rock.perm,1));
else
   error('could not open file');
end
%assert(size(rock.perm, 2)==1); % Scalar permeability
if(size(rock.perm, 2)==1)
   fprintf(fp, '%g\n', rock.perm/(milli*darcy));
elseif(any(size(rock.perm, 2)==[3,2]))
   fprintf(fp, '%g %g %g\n', rock.perm'/(milli*darcy));
else
   error('wrong dimension of rock')
end
fclose(fp);

fp = fopen(sprintf('%s-poro.dat', caseName),'w');
fprintf(fp, '%d\n', numel(rock.poro));
fprintf(fp, '%g\n', rock.poro);
fclose(fp);
end

function writeSimGrid(caseName, G)
numCells     = G.cells.num;
numCellFaces = size(G.cells.faces, 1);
numFaces     = G.faces.num;
numNodes     = G.nodes.num;
dim          = size(G.nodes.coords,2);

%% Write grid geometry specification --------------------------
%
% File format is
%   [header -- inconsequential in MATLAB]
%   [blank]
%   [numGlobalNodes]
%   [blank]
%   [numGlobalNodes lines of X Y Z coordinates]
%   [blank]
%   [numGlobalFaces]
%   [numGlobalFaces lines of unit face normals]
%   [blank]
%   [numGlobalFaces]
%   [numGlobalFaces lines of face centroids]
%   [blank]
%   [numGlobalFaces]
%   [numGlobalFaces lines of face areas]
%   [blank]
%   [numGlobalCells]
%   [numGlobalCells lines of cell centroids]
%   [blank]
%   [numGlobalCells]
%   [numGlobalCells lines of cell volumes]

filename = [caseName, '-geom.dat'];
[fid, msg] = fopen(filename, 'w');
if fid < 0, error(msg); end

%  Header
   fprintf(fid, 'geometry 0:%1d:point %1d:%1d:normal %1d:%1d:centroid %1d:1:area %1d:%1d:centroid %1d:1:volume\n', ...
		 dim, dim-1,dim, dim-1,dim, dim-1, dim,dim, dim);

%  Global nodes. XYZ-coordinates
   fprintf(fid, '\n%d\n', numNodes);
   if dim == 3,
      fprintf(fid, '%f %f %f\n', G.nodes.coords');
   else
      fprintf(fid, '%f %f\n', G.nodes.coords');
   end;

%  Unit face normals
   fprintf(fid, '\n%d\n', numFaces);
   if dim == 3,
      fprintf(fid, '%f %f %f\n', (G.faces.normals ./ G.faces.areas(:, [1, 1, 1]))');
   else
      fprintf(fid, '%f %f\n', (G.faces.normals ./ G.faces.areas(:, [1, 1]))');
   end;

%  Face centroids
   fprintf(fid, '\n%d\n', numFaces);
   if dim == 3,
      fprintf(fid, '%f %f %f\n', G.faces.centroids');
      %fprintf(fid, '%f %f %f\n', G.faces.centroids');
   else
      fprintf(fid, '%f %f\n', G.faces.centroids');
   end;

%  Face areas
   fprintf(fid, '\n%d\n', numFaces);
   fprintf(fid, '%f\n', G.faces.areas);

%  Cell centroids
   fprintf(fid, '\n%d\n', numCells);
   if dim == 3,
      fprintf(fid, '%f %f %f\n', G.cells.centroids');
   else
      fprintf(fid, '%f %f\n', G.cells.centroids');
   end;

%  Cell volumes
   fprintf(fid, '\n%d\n', numCells);
   fprintf(fid, '%f\n', G.cells.volumes);
   fclose(fid);


%-------------------------------------------------------------
%% Write grid topology specification --------------------------
%
% File format is
%   [header -- inconsequential in MATLAB]
%   [blank]
%   [numCells numCellFaces numFaces numNodes]
%   [blank]
%   [numCells lines of
%       nLocFaces  locFace(1 : nLocFaces)], locFace indexed from zero.
%   [blank]
%   [numCellFaces lines of cellFace/orientation pairs]
%   [blank]
%   [numFaces lines of
%       nNodes globNode(1 : nNodes)], globNode indexed from zero.


filename = [caseName, '-topo.dat'];
[fid, msg] = fopen(filename, 'w');
if fid < 0, error(msg); end

%  Header
fprintf(fid, 'topology %1d %1do %1d 0 %1d-%1do %1do-%1d %1d-0\n', ...
        dim, dim-1, dim-1, dim,dim-1, dim-1,dim-1, dim-1);

fprintf(fid, '\n%d %d %d %d\n\n',numCells,numCellFaces,numFaces,numNodes);


for cellIx = 1 : numCells
    range = (G.cells.facePos(cellIx) : G.cells.facePos(cellIx+1) - 1) - 1;
    fprintf(fid, '%d', size(range,2));
    fprintf(fid, ' %d', range);
    fprintf(fid, '\n');
end
fprintf(fid, '\n');

cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
sgn = -1+2*(G.faces.neighbors(G.cells.faces(:,1), 1) == cellNo);

for cellFaceIx = 1 : numCellFaces
   orientation = sgn(cellFaceIx);%G.cells.faces(cellFaceIx,3);
   if  orientation  == -1, orientation = 0; end
   fprintf(fid, '%d %d\n', G.cells.faces(cellFaceIx,1) - 1, orientation);
end
fprintf(fid, '\n');


for faceIx = 1 : numFaces
    range = (G.faces.nodePos(faceIx) : G.faces.nodePos(faceIx+1) - 1);
    fprintf(fid, '%d', size(range,2));
    fprintf(fid, ' %d', G.faces.nodes(range) - 1);
    fprintf(fid, '\n');
end

fclose(fid);
end
