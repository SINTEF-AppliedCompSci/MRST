function spill_illustration()

   [Lx,Ly,H] = deal(10000, 5000, 50);
   G  = cartGrid([100 100 1],[Lx Ly H]);
   x  = G.nodes.coords(1:G.nodes.num/2,1)/Lx; 
   y  = G.nodes.coords(1:G.nodes.num/2,2)/Ly;
   z  = G.nodes.coords(1:G.nodes.num/2,3)/H;
   zt = z + x - 0.2*sin(5*pi*x).*sin(5*pi*y.^1.5) - 0.075*sin(1.25*pi*y) + 0.15*sin(x+y);
   zb = 1 + x;
   G.nodes.coords(:,3) = [zt; zb]*H+1000;

   % making a subset of the grid, away from boundaries
   ix = reshape(1:(100^2), 100, 100);
   ix = ix(20:80, 20:80);
   G = extractSubgrid(G, ix(:));

   G = computeGeometry(G);

   Gt = topSurfaceGrid(G);
   ta  = trapAnalysis(Gt, false);
   tan = computeNodeTraps(Gt);


   %% Plot start point
   plotGrid(Gt, 'facecolor', 'white','facealpha', 0.8); hold on; axis off;
   view(76,22);
   start_node = 3082;
   plot_node(Gt, start_node); hold on;
   keyboard;
   
   %% Plot spill path
   cur_node = start_node;
   next_node = tan.dstr_neigh(cur_node);
   while next_node ~= 0
      plot_edge(Gt, cur_node, next_node);
      cur_node = next_node;
      next_node = tan.dstr_neigh(cur_node);
   end
   
   %% plot sommet point
   plot_node(Gt, cur_node);
   
   %% Plot filling trap
   
   % trap almost empty
   keyboard;
                    plot_partial_fill(G, Gt, ta, tan, cur_node, 0.2);
                    plot_partial_fill(G, Gt, ta, tan, cur_node, 0.7);
   [z1_ix, z2_ix] = plot_partial_fill(G, Gt, ta, tan, cur_node, 1);
   
   %% plot spill point
   if tan.trap_regions(z1_ix) == tan.trap_regions(cur_node)
      z_ix = z2_ix;
   else
      z_ix = z1_ix;
   end    
   plot_node(Gt, z_ix);
   keyboard;

   %% Plot second spill path
  
   cur_node = z_ix;
   next_node = tan.dstr_neigh(cur_node);
   while next_node ~= 0
      plot_edge(Gt, cur_node, next_node);
      cur_node = next_node;
      next_node = tan.dstr_neigh(cur_node);
   end

   % Plot sommet
   plot_node(Gt, cur_node);
   
   keyboard;
   %% Plot second trap filled
   plot_partial_fill(G, Gt, ta, tan, cur_node, 1);
   
end

function [spill_z1_ix, spill_z2_ix] = plot_partial_fill(G, Gt, ta, tan, sommet_node, fill_degree)

   sommet_z = Gt.nodes.z(sommet_node);
   trap_num = tan.traps(sommet_node);
   spill_edge = tan.trap_sp_edge_ix{trap_num}; 
   assert(numel(spill_edge)) = 1;
   spill_z1_ix = Gt.faces.nodes(Gt.faces.nodePos(spill_edge));
   spill_z2_ix = Gt.faces.nodes(Gt.faces.nodePos(spill_edge)+1);
   spill_z1 = Gt.nodes.z(spill_z1_ix);
   spill_z2 = Gt.nodes.z(spill_z2_ix);
   spill_z = max(spill_z1, spill_z2);
   
   
   h = zeros(Gt.cells.num, 1);
   trapcells = (ta.traps==trap_num);
   
   level_z = (1-fill_degree) * sommet_z + fill_degree * spill_z;
   
   h(trapcells) = level_z - Gt.cells.z(trapcells);
   h(h<0) = 0;
   
   hold on;
   plotPlume(G, Gt, h(:), 'facecolor', [1 .5 .5], 'facealpha', 0.8);

end


function plot_edge(Gt, n1, n2)

   c = [Gt.nodes.coords(n1, :), Gt.nodes.z(n1); ...
        Gt.nodes.coords(n2, :), Gt.nodes.z(n2)];
   plot3(c(:,1), c(:,2), c(:,3), 'r', 'linewidth', 3);
   
   
end


function plot_node(Gt, n)

   coord = [Gt.nodes.coords(n, :), Gt.nodes.z(n)];
   plot3(coord(1), coord(2), coord(3), 'ro', 'markersize', 10, ...
         'markerfacecolor', 'r');
   
end

   
   