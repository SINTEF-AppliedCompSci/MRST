function stresstest_coordmap(G)

   cmap = coordMap();
   
   max_ix = 0;
   for f = G.faces.nodes'
   
      coord = G.nodes.coords(f,:);
      
      ix = cmap(coord);
      
      max_ix = max(max_ix, ix);
      
   end
      
   max_ix
end
