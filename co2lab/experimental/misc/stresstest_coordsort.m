function stresstest_coordsort(G)


   nodes_total = zeros(numel(G.faces.nodes), 4);
   line = 1;

   % Registering all nodes, with duplicates
   tic;
   for f = G.faces.nodes'
      nodes_total(line, :) = [G.nodes.coords(f,:), line];
      line = line+1;
   end
   toc;
   
   %sorting rows
   tic; 
   nodes_total = sortrows(nodes_total);
   
   keep_ix = [true; any(diff(nodes_total(:,1:3)), 2)];
   
   nodes_total = nodes_total(keep_ix,:);
   toc;
   
   sum(keep_ix)
   size(nodes_total, 1)
   
end
