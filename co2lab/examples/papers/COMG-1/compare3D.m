function compare3D()
   
   %% Compute VE simulation result
   simres_VE = runStandardModel('data/residualExample1Data_3dcompare_VE', ...
                             @(dummy) 0                              , ...
                                'A'             , 0                     , ...
                                'residual'      , [false true]          , ...
                                'dis_types'     , {'none'}, ...
                                'xres', 200);
   
   %% Compute corresponding 3D simulation result
   simres_3D = runStandardModel3D('data/residualExample1Data_3dcompare_3D', ...
                                  @(dummy) 0                              , ...
                                  'residual'    , [false true], ...
                                  'xres', 200, ...
                                  'zres', 30);
   

   %% Produce comparison graph
   Gt = simres_VE{1}.Gt;
   xc = Gt.cells.centroids(:,1) / 1e3;
   
   tstep = numel(simres_VE{1}.states) - 70;

   % Vertical equilibrium profile, no residual saturation
   VE_nores_h = simres_VE{1}.states{tstep}.s(:,2) .* simres_VE{1}.Gt.columns.dz;

   % Vertical equilibrium profile, with residual saturation
   resvals = simres_VE{2}.fluid_params.args{end};
   VE_res_h   = free_sg(simres_VE{2}.states{tstep}.s(:,2), ...
                        simres_VE{2}.states{tstep}.sGmax, ...
                        struct('res_water', resvals(1), ...
                               'res_gas', resvals(2))) .* ...
                simres_VE{1}.Gt.columns.dz;
   
   % 3D profile, without residual saturation
   Gt2 = topSurfaceGrid(simres_3D{1}.G);
   sG = simres_3D{1}.states{tstep}.s(:,2);
   nores_h_3D = sum(reshape(sG .* col(sortrows([Gt2.columns.cells, Gt2.columns.dz], 1), 2), ...
                            numel(xc), ...
                            []), ...
                    2);
   
   % 3D profile, with residual saturation
   sG    = simres_3D{2}.states{tstep}.s(:,2);
   sGmax = simres_3D{2}.states{tstep}.sGmax;

   %h_loc = (sG - (sGmax - sG) .* resvals(2)) ./ (1 - resvals(1));
   
   % removing residual saturation
   sG((sG < sGmax) & (sG <= resvals(2))) = 0;
   res_h_3D = sum(reshape(sG .* col(sortrows([Gt2.columns.cells, Gt2.columns.dz], 1), 2), ...
                          numel(xc), ...
                          []), ...
                  2);
   
   %% Comparion plot
   figure; % new figure
   clf; hold on;
   plot(xc , VE_nores_h , 'b-'  , 'LineWidth' , 2);
   plot(xc , VE_res_h   , 'r-'  , 'LineWidth' , 2);
   plot(xc , nores_h_3D , 'b--' , 'LineWidth' , 2);
   plot(xc , res_h_3D   , 'r--' , 'LineWidth' , 2);
   
   axis tight
   set(gca, 'YDir', 'reverse', 'FontSize', 16);
   legend('VE, no residual', 'VE, residual', '3D, no residual', '3D, residual');
   
   %% Grid-based plot of 3D solution
   
   figure; % new figure
   % Straigthening grid to make it mode plot-friendly
   G = simres_3D{1}.G;
   num_layer_nodes = (G.cartDims(1) + 1) * (G.cartDims(2) + 1);
   G.nodes.coords(:,3) = ...
       G.nodes.coords(:,3) - ...
       repmat(G.nodes.coords(1:num_layer_nodes, 3), ...
              prod(G.cartDims + 1) / num_layer_nodes, 1);
       
   
   [subG, gc] = extractSubgrid(G, 1:(G.cartDims(1) * 9));

   field = simres_3D{1}.states{tstep}.s(:,2);
   
   subG.nodes.coords(:,1) = subG.nodes.coords(:,1)/1000; % in km
   plotCellData(subG, field(gc), 'edgecolor', 'none'); view(0,0); 
   colormap winter;
   set(gca, 'fontsize', 16);
   axis tight;
   zlabel('depth [m]');
   xlabel('Lateral extent [km]');
   set(gcf,'Renderer', 'painters');
end

function c = col(mat, colnum)
   c = mat(:,colnum);
end
