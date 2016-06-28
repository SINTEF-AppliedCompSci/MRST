function plotFigs3DCase(simres)
   
% Produce plot for comparison with figure 5
   Gt = topSurfaceGrid(simres{1}.G);
   xc = Gt.cells.centroids(:,1) / 1e3;
   
   % Curve without residual trapping
   %state = simres{1}.states{end-70};
   state = simres{1}.states{end-70};
   sG = state.s(:,2);
   h = sum(reshape(sG .* col(sortrows([Gt.columns.cells, Gt.columns.dz], 1), 2), ...
               numel(xc), ...
               []), ...
           2);
   plot(xc, h, 'b-', 'LineWidth', 2);
   
   axis tight;
   set(gca, 'YDir', 'reverse', 'FontSize', 16);
   legend('3D, no residual', 'location', 'east');
end

function c = col(mat, colnum)
   c = mat(:,colnum);
end
