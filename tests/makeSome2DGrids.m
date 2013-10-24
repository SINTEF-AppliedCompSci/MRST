function [g_top,g,rock,grdecl]=makeSome2DGrids(mycase)
%{
Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.

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

%mycase = 'sleipner'
switch mycase
   case {'testgrid'}
      test_grid
      %nx=40;ny=120;nz=1;
      nx=10;ny=10;nz=2;
      dims= [nx,ny,nz];
      %physdims=dims.*[3000,9000,80];
      physdims=[3000,9000,80];
      g = cartGrid(dims,physdims); %, [1 1 1]);
      %g.nodes.coords(:,[1,2])=bsxfun(@minus,g.nodes.coords(:,[1,2]),physdims(:,[1,2])/2);
      r2=(g.nodes.coords(:,1).^2+g.nodes.coords(:,2).^2);
      %g.nodes.coords(:,3)=g.nodes.coords(:,3)-10*physdims(3).*(1-r2/max(r2));
      g.nodes.coords(:,3)=2000+g.nodes.coords(:,3)-0.2*g.nodes.coords(:,1);
      %g = removeCells(g, [1 4]);x
      %g.nodes.coords = twister(g.nodes.coords);
      g = computeGeometry(g);
      %rock.perm = 100*rand(g.cells.num,1)*milli*darcy;
      rock.perm = ones(g.cells.num,1)*100*milli*darcy;
      rock.poro = 0.3*ones(g.cells.num,1);
      rock.poro(1:nx*ny) = 0.01;
      %figure
      %plotGrid(g); view(3)
      grdecl = [];
   case {'norne'}
      grdecl_file = '../ioNorne/data_io/BC0407_IO.DATA';
      grdecl = readGRDECL(grdecl_file,'verbose',true)
      GG = processGRDECL(grdecl,'Tolerance',1.0); %clear grdecl;
      g = computeGeometry(GG(1));
      rock = grdecl2Rock(grdecl,g.cells.indexMap);
      rock.perm = convertFrom(rock.perm, milli*darcy);
      %g.cells.faces = g.cellFaces;
      %g.faces.nodes = g.faceNodes;
   case {'saigup'}
      makeSAIGUP
   case {'johansen'}
      sector = 'NPD5';
      sector = fullfile(ROOTDIR, 'projects', 'co2', 'johansen', sector);
      grdecl = readGRDECL([sector, '.grdecl']);
      g = processGRDECL(grdecl); %clear grdecl;
      g = computeGeometry(g);
      K = reshape(load([sector, '_Permeability.txt'])', [], 1);
      p = reshape(load([sector, '_Porosity.txt'])',     [], 1);
      rock.perm = bsxfun(@times, [1 1 0.1], K(g.cells.indexMap)).*milli*darcy;
      rock.poro = p(g.cells.indexMap);
      clear p K;
      fluid = initSimpleFluid('mu' , [0.307 0.049] .* centi*poise, ...
                              'rho', [973 617] .* kilogram/meter^3, ...
                              'n'  , [2, 2]);
      gravity off

      %% Set well and boundary data
      % Read the well data. Injection rate is 1.4e4 m^3/day of supercritical CO2.
      w = load([sector, '_Well.txt']);
      W = verticalWell([], g, rock,  w(1,1), w(1,2), w(1,3):w(1,4),  ...
         'Type', 'rate', 'Val', 1.4e4*meter^3/day,  ...
         'Radius', 0.1, 'Comp_i', [1,0], 'name', 'I', 'comp_i', 0.0);

      %%
      % Find outer faces and set pressure boundary conditions. Here we need to do
      % some manual fiddling with the parameters to 'pside' to avoid setting
      % pressure conditions on outer boundaries that appear at faults inside
      % the lateral outline of the model.
      nx = g.cartDims(1); ny=g.cartDims(2); nz=g.cartDims(3); n = 1;
      figure;
      plotGrid(g, 'faceColor', 'none', 'EdgeAlpha', 0.1)
      ix1 = boundaryFaceIndices(g, 'BACK', 1:nx-6, 1:nz,  1:4);      
      plotFaces(g, ix1, 'r')
      ix2 = boundaryFaceIndices(g, 'BACK', nx-5:nx, 1:nz, 1:ny);
      plotFaces(g, ix2, 'g')
      ix3 = boundaryFaceIndices(g, 'RIGHT',    1:ny/2, 1:nz,  nx-4:nx);
      plotFaces(g, ix3, 'y')
      ix4 = boundaryFaceIndices(g, 'RIGHT', ny/2+1:ny, 1:nz ,  1:nx);
      ix5 = boundaryFaceIndices(g, 'FRONT', 1:nx, 1:nz, ny/2:ny);
      plotFaces(g, ix4, 'm')

      
      
      nx = g.cartDims(1); ny=g.cartDims(2); nz=g.cartDims(3); n = 1;
      bc  = pside([], g, 'LEFT',  290*barsa,      1:ny, 1:nz, 'range',    1:20, 'sat',1);
      bc  = pside(bc, g, 'RIGHT', 290*barsa,    1:ny/2, 1:nz, 'range', nx-4:nx, 'sat',1);
      bc  = pside(bc, g, 'RIGHT', 290*barsa, ny/2+1:ny, 1:nz, 'range',    1:nx, 'sat',1);
      bc  = pside(bc, g, 'BACK',  290*barsa,    1:nx-6, 1:nz, 'range',    1:4,  'sat',1);
      bc  = pside(bc, g, 'BACK',  290*barsa,   nx-5:nx, 1:nz, 'range',    1:ny, 'sat',1);
      bc  = pside(bc, g, 'FRONT', 290*barsa,      1:nx, 1:nz, 'range', ny/2:ny, 'sat',1);
      
      
      
      
      
      
      
   case {'sleipner'}
      grdecl=readGRDECL(['data/sleipner/sleipner.data'])
      %grdecl.ZCORN=-grdecl.ZCORN;
      %ind = [40,70;100,140;3,7];%a bit smaller testcase around injecto
      %ind = [60,62;125,127;3,3];%small test
      ind = [1,120;1,200;3,7];%full model without the top surface
      %sleipner has lefthanded numbering
      grdecl = cutGrdecl(grdecl,ind,'lefthanded_numbering',true);
      g=processGRDECL(grdecl);
      g=computeGeometry(g);
      rock = grdecl2Rock(grdecl,g.cells.indexMap);
      rock.perm=rock.perm*milli*darcy;
   case {'sleipner_big'}
      grdecl=readGRDECL(['/home/hnil/mywork/sleipner/grid/sleipner_refined_xy3.data'])
      %grdecl.ZCORN=-grdecl.ZCORN;
      %ind = [40,70;100,140;3,7];%a bit smaller testcase around injecto
      %ind = [60,62;125,127;3,3];%small test
      ind = [1,3*120;1,3*200;3,7];%full model without the top surface
      %sleipner has lefthanded numbering
      grdecl = cutGrdecl(grdecl,ind,'lefthanded_numbering',true);
      g=processGRDECL(grdecl);
      g=computeGeometry(g);
      rock = grdecl2Rock(grdecl,g.cells.indexMap);
      rock.perm=rock.perm*milli*darcy;
   otherwise
      error('No such grid')
end
%%
%figure(1),clf,plotGrid(g,'FaceColor','none','EdgeAlpha',0.2)
%figure(2),clf,plotCellData(g,rock.perm(:,1));
%%
%figure(3)
g_top=topSurfaceGrid(g, 'grdecl', grdecl);
%plotGrid(g_top)
%plotCellData(g_top, g_top.cells.z)
%plotCellData(g, g.cells.centroids(:,3))
