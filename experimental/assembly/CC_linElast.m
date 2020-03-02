function [uu, out] = CC_linElast(G, C, el_bc, load, varargin)
% solve linear elastisity problem using the cell centered Multi Point Stress
% approximation by Jan Nodbotten. This code wrapps his original code.
% It has for now limited support for boundary conditions
%
% SYNOPSIS:
%   [uu, S, A, extra] = VEM2D_linElast(G, E, nu, el_bc, load, varargin)
%
% OUTPUT
%    uu is matrix of size [G.nodes.num, G.griddim] of displace ment on
%        nodes
%
% PARAMETERS:
%   G       - Grid structure as described by grid_structure, which all so
%             has the mappings from G = mrstGridWithFullMappings(G). Some options
%             need  G = computeGeometryCalc(G).
%   E       - Youngs modulo  %http://en.wikipedia.org/wiki/Shear_modulus
%   nu      - Poisson's ratio
%   el_bc   -  boundary condition of type struct('disp_bc', [], 'force_bc', [])
%
%
%
%
%   'pn'/pv - List of 'key'/value pairs for supplying optional parameters.
%             The supported options are
%            opt = struct('methodS', 'U');, ... % method 
%            
%
%
% OPTIONAL PARAMETERS:
%
%
%
% COMMENTS:
%   For all the solvetypes exept lorenzo the assembly code is for general 
%   elastisty tensor in voit's notation http://en.wikipedia.org/wiki/Voigt_notation
% SEE ALSO:
% square_2D_example, convergens_test, plotNodeData, plotGridDeformed


%{ 
Copyright 2014 Jan Nordbotten jan.nordbotten@math.uib.no
%} 
%opt = struct('methodS', 'AU');
%opt = merge_options(opt, varargin{:});
%  convert to Jan Norbottens grid structure based on cell arrays
   opt=struct('pressure',[]);
   opt=merge_options(opt,varargin{:});
   gem = mrstGrid2NG(G);

   BCMass = 'Dirichlet';
   BCMoment = 'Dirichlet';
   isbndMass = gem.isbnd;
   isbndMoment = isbndMass;
   [isbndMoment{strcmp({isbndMass{:}}, BCMass)}] = deal('Neumann');

   for iter1 = find(~strcmp({gem.isbnd{:}}, 'no'))
      isbndMoment{iter1} = BCMoment;
   end

   gem.face_tangent = []; % not used maybe only for definition of stress??
                          % Stress discretization
   gem.isbnd = isbndMoment;
   isdirec = zeros(G.faces.num, 1);
   isdirec(el_bc.disp_bc.faces) = 1;
   isneum = zeros(G.faces.num, 1);
   if(~all(el_bc.disp_bc.mask(:)>0))
      warning('Partial Direchlet boundary conditions not supported for now')
   end


   if(~isempty(el_bc.force_bc))
      isneum(el_bc.force_bc.faces) = 1;
   end

   % correct boundary condition type before assembly
   % NB there is some simplification approximation for direclet boudnary
   %   so simple patch test do not work exactly
   for i = 1:G.faces.num
      if(all(G.faces.neighbors(i, :)>0))
         assert(strcmp(gem.isbnd{i}, 'no'))
      else
         if(isdirec(i))
            assert(~isneum(i))
            gem.isbnd{i} = 'Dirichlet';
         else
            gem.isbnd{i} = 'Neumann'; 
         end
         
      end
   end

   % define rock properties
   constit = voigtToFull(G, C);
   %out = MPSA(gem, constit, 'FV', 1);
   gem.G=G;
   %out = MPSA(gem, constit, 'FV', 1);
   out = MPSA_stab(gem, constit);
   
   %md = MPSA_vectorized(G,constit,'weakCont',weakCont,'bc',bc,'eta',eta,'biot',0);
   %%{
   %isBoundary = any(G.faces.neighbors == 0,2);
   %bc = addBC([],find(isBoundary),'pressure',0);
   %out = MPSA_vectorized(G,constit,'weakCont',1,'bc',bc,'eta',0,'biot',0);
   %out.rhsA=1;
   %}
   out.constit=constit;
   
   % 
   %moment = out.A;
   %stress = out.stress;
   %divD = out.divD;
   %rhsA = out.rhsA;


   bc = el_bc.disp_bc;
   bcfaces = bc.faces;
   val = bc.uu_face;
   bc_disp = zeros(G.faces.num, G.griddim);
   bc_cond = zeros(G.faces.num, 1);
   bc_cond(bcfaces, :) = 1;
   bc_disp(bcfaces, :) = val;
   cellno = rldecode([1:G.cells.num]', diff(G.cells.facePos));%#ok
                                                              % calculate mean of displace ment of bounday conditions
   bc_cell_disp = zeros(G.cells.num, G.griddim);
   for i = 1:G.griddim
      bc_cell_disp(:, i) = accumarray(cellno, bc_disp(G.cells.faces(:, 1), i));
   end
   nbc = accumarray(cellno, bc_cond(G.cells.faces(:, 1)));
   bc_cell_disp(nbc>0, :) = bsxfun(@rdivide, bc_cell_disp(nbc>0, :), nbc(nbc>0));

   bc_cell_disp = reshape(bc_cell_disp', [], 1);
   bc_disp = reshape(bc_cell_disp', [], 1);

   % forces on each cell
   force = -bsxfun(@times, load(G.cells.centroids), G.cells.volumes);

   %
   if(~isempty(el_bc.force_bc))
      force_bc = el_bc.force_bc;
      bccell = sum(G.faces.neighbors(force_bc.faces, :), 2);
      area = G.faces.areas(force_bc.faces);
      bcforce_f = bsxfun(@times, force_bc.force, area);
      % avoid problem repititions
      bcforce = zeros(G.cells.num, G.griddim);
      for i = 1:G.griddim
         bcforce(:, i) = accumarray(bccell, bcforce_f(:, i), [G.cells.num, 1]);
      end
      force = force + bcforce;
   end
   force = force';
   force = force(:);
   % The Diriclet boundary conditions momentRHS.A*bc_disp) is cell based and
   % is not accurate if one cell has more than 1 boundary neighbor
   % rhs = force+momentRHS.A*bc_dips

   rhs = force + out.rhsA*bc_disp;
   if(~isempty(opt.pressure))
       rhs=rhs+out.gradP*opt.pressure;
   end
   
   uu = out.A\rhs;
   %uu = moment\(force+(rhsA*bc_disp));%
   uu = reshape(uu, G.griddim, G.cells.num)';

end