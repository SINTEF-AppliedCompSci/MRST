function table = reproduceDensityTable()

   CO2 = CO2props('rho_huge','');
   
   orig = CO2.rhoInfo();
   figure(1); surf(orig.vals,'edgecolor','none'); view(0, 90); 

   p_points = linspace(orig.P.span(1), orig.P.span(2), orig.P.num);
   t_points = linspace(orig.T.span(1), orig.T.span(2), orig.T.num);
   
   [tmesh, pmesh] = meshgrid(t_points, p_points);
   
   new_vals = reshape(CO2.rho(pmesh(:), tmesh(:)), orig.P.num, orig.T.num);
   figure(2); surf(new_vals, 'edgecolor','none'); view(0, 90); 
   
   figure(3);surf(new_vals - orig.vals, 'edgecolor','none'); view(0, 90);
   
end