function G = squareGrid(cartDims, L, varargin)
% G = squareGrid(cartDims, L, varargin)
% make square with different grid types starting from cartesian
% OPTIONS

%opt = struct('L', [1 1], ...
%    'cartDims', [20 20], ...
%    'grid_type', 'square', ...
%    'disturb', 0.0);
   opt = struct('grid_type', 'cartgrid', ...
                'disturb', 0.0);
   opt = merge_options(opt, varargin{:});
   opt.L = L;opt.cartDims = cartDims;
   if(numel(opt.cartDims) == 2)
      [nx, ny] = deal(opt.cartDims(1), opt.cartDims(2));[Lx, Ly] = deal(opt.L(1), opt.L(2));
   else
      assert(strcmp(opt.grid_type, 'cartgrid'))
   end
   switch opt.grid_type
     case 'cartgrid'
       G = cartGrid(opt.cartDims, opt.L);
       if(opt.disturb>0)
          G = twister(G, opt.disturb);
       end
       %G = computeGeometry(G);
     case 'triangle'
       G_tmp = cartGrid(opt.cartDims, opt.L);
       if(opt.disturb>0)
          G_tmp = twister(G_tmp, opt.disturb);
       end
       p = [G_tmp.nodes.coords(:, 1), G_tmp.nodes.coords(:, 2)];
       t = delaunayn(p);
       G = triangleGrid(p, t);
       %G = computeGeometry(G);
     case 'triangle2'  
        rng('default')
       fd=@(p) drectangle(p,0,opt.L(1),0,opt.L(2));
       fh=@(p) huniform(p);
       [p,t]=distmesh2d(fd,fh,sqrt(prod(opt.L)/prod(opt.cartDims)),[0,0;opt.L(1),opt.L(2)],[0,0;opt.L(1),0;0,opt.L(2);opt.L(1),opt.L(2)]);
       G = triangleGrid(p, t);
       if(opt.disturb>0)
           G=twister(G,opt.disturb);
       end
      case 'triangle3'  
        rng('default')
       fd=@(p) drectangle(p,0,1,0,1);
       fh=@(p) huniform(p);
       [p,t]=distmesh2d(fd,fh,sqrt(1/prod(opt.cartDims)),[0,0;1,1],[0,0;1,0;0,1;1,1]);
       G = triangleGrid(p, t);
       if(opt.disturb>0)
           G=twister(G,opt.disturb);
       end
      G.nodes.coords=bsxfun(@times,G.nodes.coords,opt.L);  
     case 'pebi'
       G_tmp = cartGrid(opt.cartDims, opt.L);
       if(opt.disturb>0)
          G_tmp = twister(G_tmp, opt.disturb);
       end
       p = [G_tmp.nodes.coords(:, 1), G_tmp.nodes.coords(:, 2)];
       t = delaunayn(p);
       G = triangleGrid(p, t);
       G = pebi(G);
       G = removeShortEdges(G, min(opt.L)/(max(opt.cartDims)*10));
       G = sortEdges(G);
     case 'pebi2'
       rng('default')
       G_tmp = cartGrid(opt.cartDims, opt.L);
       internal = all((G_tmp.nodes.coords>0) & (G_tmp.nodes.coords< repmat(opt.L, G_tmp.nodes.num, 1)), 2);
       G_tmp.nodes.coords(internal, :) = G_tmp.nodes.coords(internal, :)...
           +bsxfun(@times, rand(size(G_tmp.nodes.coords(internal, :))), opt.disturb*opt.L./max(opt.cartDims));
       if(opt.disturb>0)
          G_tmp = twister(G_tmp, opt.disturb);
       end
       p = [G_tmp.nodes.coords(:, 1), G_tmp.nodes.coords(:, 2)];
       t = delaunayn(p);
       G = triangleGrid(p, t);
       G = pebi(G);
       G = removeShortEdges(G, min(opt.L)/(max(opt.cartDims)*10));
       G = sortEdges(G);
     case 'pebi3'
       rng('default')
       fd=@(p) drectangle(p,0,opt.L(1),0,opt.L(2));
       fh=@(p) huniform(p);
       [p,t]=distmesh2d(fd,fh,sqrt(prod(opt.L)/prod(opt.cartDims)),[0,0;opt.L(1),opt.L(2)],[0,0;opt.L(1),0;0,opt.L(2);opt.L(1),opt.L(2)]);
       G = triangleGrid(p, t);
       G = pebi(G);
       if(opt.disturb>0)
           G=twister(G,opt.disturb);
       end
       G = removeShortEdges(G, min(opt.L)/(max(opt.cartDims)*10));
       G = sortEdges(G);   
     case 'pebi4'
       rng('default')
       fd=@(p) drectangle(p,0,1,0,1);
       fh=@(p) huniform(p);
       [p,t]=distmesh2d(fd,fh,sqrt(1/prod(opt.cartDims)),[0,0;1,1],[0,0;1,0;0,1;1,1])
       [p,t,pix]=fixmesh(p,t);
       G = triangleGrid(p, t);
       G=sortEdges(G);
       G = pebi(G);
       if(opt.disturb>0)
           G=twister(G,opt.disturb);
       end
       G.nodes.coords=bsxfun(@times,G.nodes.coords,opt.L);
       G = removeShortEdges(G, min(opt.L)/(max(opt.cartDims)*10));
       G = sortEdges(G);    
     case   {'boxes1', 'boxes2', 'boxes3', 'boxes4'}
       G1 = cartGrid(opt.cartDims, opt.L);
       G2 = cartGrid(2*opt.cartDims+1, opt.L);
       G = glue2DGrid(G1, translateGrid(G2, [opt.L(1), -0.0]));
       %if strcmp(opt.grid_type, 'mixed')
       if(str2num(opt.grid_type(end))>1)
          G3 = cartGrid(opt.cartDims+2, opt.L);
          %G = glue2DGrid(G1, translateGrid(G2, [opt.L(1), -0.0]));
          G = glue2DGrid(G, translateGrid(G3, 2*[opt.L(1), 0.0]));
          if(str2num(opt.grid_type(end))>2)
             G4 = cartGrid(opt.cartDims+[7, 1], [3*opt.L(1), opt.L(2)]);
             dx = 0.05;
             G5 = cartGrid([opt.cartDims(1)*5, 1], [3*opt.L(1), dx]);
             G = glue2DGrid(G, translateGrid(G5, [0, opt.L(2)]));
             G = glue2DGrid(G, translateGrid(G4, [0, opt.L(2)+dx]));
             G = glue2DGrid(G, translateGrid(G4, [0, -opt.L(2)]));
             if(str2num(opt.grid_type(end))>3)
                G6 = cartGrid([3, 10], [opt.L(1), 3*opt.L(2)+dx]);
                %G = G1;
                %G = glue2DGrid(G1, translateGrid(G2, [opt.L(1), -0.0]));
                %G = glue2DGrid(G2, translateGrid(G1, [opt.L(1), -0.0]));
                %G = glue2DGrid(G, translateGrid(G3, 2*[opt.L(1), 0.0]));
                %G = glue2DGrid(G, translateGrid(G5, [0, opt.L(2)]));
                %G = glue2DGrid(G, translateGrid(G4, [0, opt.L(2)+dx]));
                %G = glue2DGrid(G, translateGrid(G4, [0, -opt.L(2)]));
                G = glue2DGrid(G, translateGrid(G6, [- opt.L(1), -opt.L(2)]));
                G = glue2DGrid(G, translateGrid(G6, [3*opt.L(1), -opt.L(2)]));
                %G.cells = rmfield(G.cells, 'indexMap');
             end
          end
       end
       G.cells = rmfield(G.cells, 'indexMap');
       %G = glue2DGrid(G, G4);%translateGrid(G4, [0, 0.0*opt.L(2)]));
       if(opt.disturb>0)
          G = twister(G, opt.disturb);
       end
       G.nodes.coords(:, 1) = G.nodes.coords(:, 1)-min(G.nodes.coords(:, 1));
       G.nodes.coords(:, 2) = G.nodes.coords(:, 2)-min(G.nodes.coords(:, 2));
       %G.nodes.coords(:, 2) = G.nodes.coords(:, 2)*opt.L(2)./Ly;
       Lx = max(G.nodes.coords(:, 1));
       Ly = max(G.nodes.coords(:, 2));
       G.nodes.coords(:, 1) = G.nodes.coords(:, 1)*opt.L(1)./Lx;
       G.nodes.coords(:, 2) = G.nodes.coords(:, 2)*opt.L(2)./Ly;
       
       %clf, plotGrid(G);%, plotGrid(G4, 'FaceColor', 'r')
       G = sortEdges(G);
       %error('There is repeted nodes')
     case   {'mixed1', 'mixed2', 'mixed3', 'mixed4'}
       G1 = squareGrid(opt.cartDims, opt.L, 'grid_type', 'cartgrid', 'disturb', opt.disturb);
       %G2 = cartGrid(2*opt.cartDims+1, opt.L);
       G2 = squareGrid(opt.cartDims+1, opt.L, 'grid_type', 'pebi', 'disturb', opt.disturb);
       G = glue2DGrid(G1, translateGrid(G2, [opt.L(1), -0.0]));
       %if strcmp(opt.grid_type, 'mixed')
       if(str2num(opt.grid_type(end))>1)
          G3 = squareGrid(opt.cartDims+2, opt.L, 'grid_type', 'triangle', 'disturb', 0.02);
          %G = glue2DGrid(G1, translateGrid(G2, [opt.L(1), -0.0]));
          G = glue2DGrid(G, translateGrid(G3, 2*[opt.L(1), 0.0]));
          if(str2num(opt.grid_type(end))>2)
             G4 = cartGrid(opt.cartDims+[7, 1], [3*opt.L(1), opt.L(2)]);
             %G4 = squareGrid(opt.cartDims+[7, 1], [3*opt.L(1), opt.L(2)],'grid_type', 'triangle3', 'disturb', 0.02);
             %G4 = squareGrid(opt.cartDims+[7, 1], [3*opt.L(1), opt.L(2)], 'grid_type', 'triangle', 'disturb', opt.disturb)
             % using G44 as pebi or triangel instead of G4 fails
             %G44 = squareGrid(opt.cartDims+[7, 1], [3*opt.L(1), opt.L(2)], 'grid_type', 'pebi', 'disturb', opt.disturb)
             dx = 0.05;
             G5 = cartGrid([opt.cartDims(1)*5, 1], [3*opt.L(1), dx]);
             G = glue2DGrid(G, translateGrid(G5, [0, opt.L(2)]));
             G = glue2DGrid(G, translateGrid(G4, [0, opt.L(2)+dx]));
             G = glue2DGrid(G, translateGrid(G4, [0, -opt.L(2)]));
             if(str2num(opt.grid_type(end))>3)
                %G7 = cartGrid([3, 10], [opt.L(1), 3*opt.L(2)+dx]);
                G7 = squareGrid([3, 10], [opt.L(1), 3*opt.L(2)+dx], 'grid_type', 'triangle', 'disturb', opt.disturb);
                %pebi fails if G44 is used
                G6 = squareGrid([3, 10], [opt.L(1), 3*opt.L(2)+dx], 'grid_type', 'pebi', 'disturb', opt.disturb);
                %G = G1;
                %G = glue2DGrid(G1, translateGrid(G2, [opt.L(1), -0.0]));
                %G = glue2DGrid(G2, translateGrid(G1, [opt.L(1), -0.0]));
                %G = glue2DGrid(G, translateGrid(G3, 2*[opt.L(1), 0.0]));
                %G = glue2DGrid(G, translateGrid(G5, [0, opt.L(2)]));
                %G = glue2DGrid(G, translateGrid(G4, [0, opt.L(2)+dx]));
                %G = glue2DGrid(G, translateGrid(G4, [0, -opt.L(2)]));
                G = glue2DGrid(G, translateGrid(G6, [- opt.L(1), -opt.L(2)]));
                G = glue2DGrid(G, translateGrid(G7, [3*opt.L(1), -opt.L(2)]));
                %G.cells = rmfield(G.cells, 'indexMap');
             end
          end
       end
       G.cells = rmfield(G.cells, 'indexMap');
       %G = glue2DGrid(G, G4);%translateGrid(G4, [0, 0.0*opt.L(2)]));
       if(opt.disturb>0)
          G = twister(G, opt.disturb);
       end
       G.nodes.coords(:, 1) = G.nodes.coords(:, 1)-min(G.nodes.coords(:, 1));
       G.nodes.coords(:, 2) = G.nodes.coords(:, 2)-min(G.nodes.coords(:, 2));
       %G.nodes.coords(:, 2) = G.nodes.coords(:, 2)*opt.L(2)./Ly;
       Lx = max(G.nodes.coords(:, 1));
       Ly = max(G.nodes.coords(:, 2));
       G.nodes.coords(:, 1) = G.nodes.coords(:, 1)*opt.L(1)./Lx;
       G.nodes.coords(:, 2) = G.nodes.coords(:, 2)*opt.L(2)./Ly;
       
       %clf, plotGrid(G);%, plotGrid(G4, 'FaceColor', 'r')
       G = sortEdges(G);
       %error('There is repeted nodes')
       
       
     otherwise
       error('no such grid type')
   end
   if(~strcmp(opt.grid_type, 'cartgrid'))
      Lx = max(G.nodes.coords(:, 1));
      Ly = max(G.nodes.coords(:, 2));
      tag = zeros(G.faces.num, 1);
      tol = 1e-3/max(opt.cartDims)*min(opt.L);
      G = computeGeometry(G);
      tag(abs(G.faces.centroids(:, 1))<tol) = 1;
      tag(abs(G.faces.centroids(:, 1)-Lx)<tol) = 2;
      tag(abs(G.faces.centroids(:, 2))<tol) = 3;
      tag(abs(G.faces.centroids(:, 2)-Ly)<tol) = 4;
      if(G.griddim == 3)
         Lz = max(G.nodes.coords(:, 3));
         tag(abs(G.faces.centroids(:, 3))<tol) = 3;
         tag(abs(G.faces.centroids(:, 3)-Ly)<tol) = 4;
      end
      G.cells.faces = [G.cells.faces(:, 1), tag(G.cells.faces(:, 1))];
   end

   G = mrstGridWithFullMappings(G);
   G = computeGeometryCalc(G);
end