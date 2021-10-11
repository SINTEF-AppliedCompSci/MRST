function [G, bc, test_cases] = squareTest(varargin)
% Different test cases for linear elasticity on square domains
%
% SYNOPSIS:
%   function [G, bc, test_cases] = squareTest(varargin)
%
% DESCRIPTION: Set up different test cases to test the linear elasticity
% code.
%
% OPTIONAL PARAMETERS:
%     'L'              - Physical dimension
%     'cartDims'       - Cartesian dimension
%     'grid_type'      - grid type (see squareGrid function)
%     'disturb'        - Disturbance coefficient
%     'E'              - Young's modulus
%     'nu'             - Poisson's ratio
%     'make_testcases' - If true, compute linear displacement test cases.
%     'test_name'      - test case name. The alternatives are
%                           'original_2D_test'
%                           'hanging_rod_2D'                     
%                           'weighted_slop_rod_2D_dir'           
%                           'weighted_slop_rod_2D'               
%                           'weighted_slop_rod_2D_nobcleftright' 
%                           'rod_2D_pressure'                    
%                           'rod_2D_pressure_nobcleftright'      
%                           'rod_2D_point_pressure'              
%                           'rod_2D_point_pressure_noleftrighbc' 
%                           'all'
%
% RETURNS:
%   G          - Grid structure
%   bc         - Boundary conditions
%   test_cases - Cell structure with the test cases

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
   
   opt = struct('L', [1 1], ...
                'cartDims', [20 20], ...
                'grid_type', 'cartgrid', ...
                'disturb', 0.0, ...
                'E', 1, ...  %youngs modolo
                'nu', 0.3, ...% poiso ratio
                'make_testcases', true, ...
                'test_name', 'all',...
                'triangulate',true);
   opt = merge_options(opt, varargin{:});

   if(numel(opt.cartDims) == 2)
      [nx, ny] = deal(opt.cartDims(1), opt.cartDims(2));[Lx, Ly] = deal(opt.L(1), opt.L(2));
   else
      assert(strcmp(opt.grid_type, 'cartgrid'))
   end

   G = squareGrid(opt.cartDims, opt.L, 'grid_type', opt.grid_type, 'disturb', opt.disturb);
   if(opt.triangulate)
      faces=[1:G.faces.num]';
      G=triangulateFaces(G,faces);
   end
   G = createAugmentedGrid(G);
   G = computeGeometry(G);
   bc = cell(2*G.griddim, 1);
   
   if (strcmp(opt.grid_type, 'cartgrid'))
      if(G.griddim == 2)
         oside = {'Left', 'Right', 'Back', 'Front'};       
      else
         oside = {'Left', 'Right', 'Back', 'Front', 'Bottom', 'Top'};
      end
      for i = 1:numel(oside);
         bc{i} = pside([], G, oside{i}, 0);
         bc{i} = rmfield(bc{i}, 'type');
         bc{i} = rmfield(bc{i}, 'sat');
      end
   else
      assert(G.griddim == 2);
      x = [0, Lx];
      for i = 1:2
         faces = find(abs(G.faces.centroids(:, 1)-x(i))<sqrt(eps));
         bc{i} = addBC([], faces, 'pressure', 0);
         bc{i} = rmfield(bc{i}, 'type');
         bc{i} = rmfield(bc{i}, 'sat');
      end
      y = [0, Ly];
      for i = 1:2
         faces = find(abs(G.faces.centroids(:, 2)-y(i))<sqrt(eps));
         bc{i+2} = addBC([], faces, 'pressure', 0);
         bc{i+2} = rmfield(bc{i+2}, 'type');
         bc{i+2} = rmfield(bc{i+2}, 'sat');
      end
   end

   for i = 1:numel(bc);
      inodes = mcolon(G.faces.nodePos(bc{i}.face), G.faces.nodePos(bc{i}.face+1)-1);
      nodes = unique(G.faces.nodes(inodes));
      bc{i}.el_bc = struct('disp_bc', struct('nodes', nodes, 'uu', 0, 'faces', bc{i}.face, 'uu_face', 0, 'mask', true(numel(nodes), G.griddim)), ...
                           'force_bc', []);
   end
   %
   test_cases = {};
   Ev = repmat(opt.E, G.cells.num, 1); 
   nuv = repmat(opt.nu, G.cells.num, 1);
   C = Enu2C(Ev, nuv, G);

   if(opt.make_testcases)
      % make cell array of function of linear displacement
      [sol, names] = linearDisplacement(G.griddim);
      load = @(cc) zeros(size(cc));
      
      for i = 1:numel(sol)
         test_cases{i} = makeDirTestCase(G, bc, C, load, sol{i}, sol{i}); %#ok
         test_cases{i}.name = names{i}; %#ok
      end
      
   end

   if(G.griddim == 2)
      
      switch opt.test_name
         
        case {'original_2D_test', 'all'}
          
          load = @(x) repmat([0, 1], size(x, 1), 1);
          bcdisp = @(x) [sin(x(:, 1)*pi)*0.1, sin(x(:, 1)*pi)*0.5];
          bcdisp = @(x) 0*x;
          solution = [];
          test_cases{end+1} = makeDirTestCase(G, bc, C, load, bcdisp, solution);
          test_cases{end}.name = 'original_2D_test';
          
        case {'hanging_rod_2D', 'all'}
          
          % hanging rod case
          load = @(x) repmat([0, 1], size(x, 1), 1);
          bcdisp = @(x) [sin(x(:, 1)*pi)*0.1, sin(x(:, 1)*pi)*0.5]*0.0;
          solution = [];
          bc_left_right = {bc{1:2}};%#ok
          test_cases{end+1} = makeDirTestCase(G, bc_left_right, C, load, bcdisp, solution);
          test_cases{end}.name = 'hanging_rod_2D';
    
          load=@(x) repmat([0,1],size(x,1),1)*0.1;
          bcdisp=@(x) [sin(x(:,1)*pi)*0.1,sin(x(:,1)*pi)*0.5]*0.0;
          solution=[];
          bc_left_right={bc{1}};%#ok
          test_cases{end+1}=makeDirTestCase(G,bc_left_right,C,load,bcdisp,solution);
          test_cases{end}.name='hanging_cantilever_2D';
          
        case {'weighted_slop_rod_2D_dir', 'all'}

          load = @(x) -repmat([0, 1], size(x, 1), 1);
          bcdisp = @(x) [sin(x(:, 1)*pi)*0.1, sin(x(:, 1)*pi)*0.5]*0.0;
          solution = [];
          bc_left_right = {bc{[1, 2, 3]}};%#ok
          test_cases{end+1} = makeDirTestCase(G, bc_left_right, C, load, bcdisp, solution);
          test_cases{end}.name = 'weighted_slop_rod_2D_dir';
          
        case {'weighted_slop_rod_2D', 'all'}
          
          load = @(x) -repmat([0, 1], size(x, 1), 1);
          bcdisp = @(x) [sin(x(:, 1)*pi)*0.1, sin(x(:, 1)*pi)*0.5]*0.0;
          solution = [];
          clear bc_left_right
          bc_left_right{1} = bc{1};
          bc_left_right{1}.el_bc.disp_bc.mask(:, 2) = false;
          bc_left_right{2} = bc{2};
          bc_left_right{2}.el_bc.disp_bc.mask(:, 2) = false;
          bc_left_right{3} = bc{3};
          bc_left_right{3}.el_bc.disp_bc.mask(:, 1) = false;
          test_cases{end+1} = makeDirTestCase(G, bc_left_right, C, load, bcdisp, solution);
          test_cases{end}.name = 'weighted_slop_rod_2D';
          
        
        case {'weighted_slop_rod_2D_nobcleftright', 'all'}
          
          load = @(x) -repmat([0, 1], size(x, 1), 1);
          bcdisp = @(x) [sin(x(:, 1)*pi)*0.1, sin(x(:, 1)*pi)*0.5]*0.0;
          solution = [];
          clear bc_left_right
          bc_left_right{1} = bc{3};
          bc_left_right{1}.el_bc.disp_bc.mask(:, 1) = false;
          test_cases{end+1} = makeDirTestCase(G, bc_left_right, C, load, bcdisp, solution);
          test_cases{end}.name = 'weighted_slop_rod_2D_nobcleftright';
          
          
        case {'rod_2D_pressure', 'all'}
          
          load = @(x) -repmat([0, 0], size(x, 1), 1);
          bcdisp = @(x) [sin(x(:, 1)*pi)*0.1, sin(x(:, 1)*pi)*0.5]*0.0;
          solution = [];
          clear bc_left_right
          bc_left_right{1} = bc{1};
          bc_left_right{1}.el_bc.disp_bc.mask(:, 2) = false;
          bc_left_right{2} = bc{2};
          bc_left_right{2}.el_bc.disp_bc.mask(:, 2) = false;
          bc_left_right{3} = bc{3};
          bc_left_right{3}.el_bc.disp_bc.mask(:, 1) = false;
          force = 1;
          face_force = @(x) force*ones(size(x(:, 1)));
          faces = bc{4}.face;
          force_bc = struct('faces', faces, 'force', bsxfun(@times, face_force(G.faces.centroids(faces, :)), [0 -1]));
          test_cases{end+1} = makeDirTestCase(G, bc_left_right, C, load, bcdisp, solution);
          test_cases{end}.el_bc.force_bc = force_bc;
          test_cases{end}.name = 'rod_2D_pressure';
          
          
        case {'rod_2D_pressure_nobcleftright', 'all'}
          
          load = @(x) -repmat([0, 0], size(x, 1), 1);
          bcdisp = @(x) [sin(x(:, 1)*pi)*0.1, sin(x(:, 1)*pi)*0.5]*0.0;
          solution = [];
          clear bc_left_right
          bc_left_right{1} = bc{3};
          bc_left_right{1}.el_bc.disp_bc.mask(:, 1) = false;
          force = 1; 
          face_force = @(x) force*ones(size(x(:, 1)));
          faces = bc{4}.face;
          force_bc = struct('faces', faces, 'force', bsxfun(@times, face_force(G.faces.centroids(faces, :)), [0 -1]));
          test_cases{end+1} = makeDirTestCase(G, bc_left_right, C, load, bcdisp, solution);
          test_cases{end}.el_bc.force_bc = force_bc;
          test_cases{end}.name = 'rod_2D_pressure_nobcleftright';
          
        case {'rod_2D_point_pressure', 'all'}
          
          load = @(x) -repmat([0, 0], size(x, 1), 1);
          bcdisp = @(x) [sin(x(:, 1)*pi)*0.1, sin(x(:, 1)*pi)*0.5]*0.0;
          solution = [];
          clear bc_left_right
          bc_left_right{1} = bc{1};
          bc_left_right{1}.el_bc.disp_bc.mask(:, 2) = false;
          bc_left_right{2} = bc{2};
          bc_left_right{2}.el_bc.disp_bc.mask(:, 2) = false;
          bc_left_right{3} = bc{3};
          bc_left_right{3}.el_bc.disp_bc.mask(:, 1) = false;
          sigma = Ly/10;
          force = 1e-2;
          face_force = @(x) force*exp(-(((x(:, 1)-Lx/2))./sigma).^2);
          faces = bc{4}.face;
          force_bc = struct('faces', faces, 'force', bsxfun(@times, face_force(G.faces.centroids(faces, :)), [0 -1]));
          test_cases{end+1} = makeDirTestCase(G, bc_left_right, C, load, bcdisp, solution);
          test_cases{end}.el_bc.force_bc = force_bc;
          test_cases{end}.name = 'rod_2D_point_pressure';
          
        case {'rod_2D_point_pressure_noleftrighbc', 'all'}
          
          load = @(x) -repmat([0, 0], size(x, 1), 1);
          bcdisp = @(x) [sin(x(:, 1)*pi)*0.1, sin(x(:, 1)*pi)*0.5]*0.0;
          solution = [];
          clear bc_left_right
          bc_left_right{1} = bc{3};
          bc_left_right{1}.el_bc.disp_bc.mask(:, 1) = false;
          sigma = Ly/10;
          force = 1;
          face_force = @(x) force*exp(-(((x(:, 1)-Lx/2))./sigma).^2);
          faces = bc{4}.face;
          force_bc = struct('faces', faces, 'force', bsxfun(@times, face_force(G.faces.centroids(faces, :)), [0 -1]));
          test_cases{end+1} = makeDirTestCase(G, bc_left_right, C, load, bcdisp, solution);
          test_cases{end}.el_bc.force_bc = force_bc;
          test_cases{end}.name = 'rod_2D_point_pressure_noleftrighbc';
          
        otherwise

          error('test name not recognized.');
      end
      
   else
      assert(G.griddim == 3)
      load = @(x) repmat([0, 0, 1], size(x, 1), 1);
      bcdisp = @(x) zeros(size(x));
      solution = [];
      test_cases{end+1} = makeDirTestCase(G, bc, C, load, bcdisp, solution);
      test_cases{end}.name = 'original_3D_test';
      
      load = @(x) repmat([0, 0, 1], size(x, 1), 1);
      bcdisp = @(x) zeros(size(x));
      solution = [];
      bc_left_right = {bc{1:4}};%#ok
      test_cases{end+1} = makeDirTestCase(G, bc_left_right, C, load, bcdisp, solution);
      test_cases{end}.name = 'hanging_plate_3D';
   end


end

function testcase = makeDirTestCase(G, bc, C, load, bcdisp, solution)
   nodes = [];
   faces = [];
   mask = [];
   for i = 1:numel(bc)
      nodes = [nodes; bc{i}.el_bc.disp_bc.nodes];%#ok
      faces = [faces; bc{i}.el_bc.disp_bc.faces];%#ok
      mask = [mask; bc{i}.el_bc.disp_bc.mask];%#ok
   end
   disp_node = bcdisp(G.nodes.coords(nodes, :));
   disp_faces = bcdisp(G.faces.centroids(faces, :));
   el_bc = struct('disp_bc', struct('nodes', nodes, 'uu', disp_node, 'faces', faces, 'uu_face', disp_faces, 'mask', mask), ...
                  'force_bc', []);
   testcase = struct('el_bc', el_bc, 'C', C, 'load', load, 'disp_sol', bcdisp, 'solution', solution);
end


