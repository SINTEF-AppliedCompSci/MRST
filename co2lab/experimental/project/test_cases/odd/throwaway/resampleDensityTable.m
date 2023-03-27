function res = resampleDensityTable(pspan, tspan, table, num_psamples, num_tsamples)

   CO2 = CO2props(table,'');
   
   p_points = linspace(pspan(1), pspan(2), num_psamples);
   t_points = linspace(tspan(1), tspan(2), num_tsamples);

   [tgrid, pgrid] = meshgrid(t_points, p_points);
   
   % computing new table of values
   vals = reshape(CO2.rho(pgrid(:), tgrid(:)), num_psamples, num_tsamples);

   % setting information about P 
   P.span = pspan;
   P.num = num_psamples;
   P.stepsize = diff(P.span) / (P.num - 1);
   
   % setting information about T
   T.span = tspan;
   T.num = num_tsamples;
   T.stepsize = diff(T.span) / (T.num - 1);

   % saving result
   save rho_upsampled P T vals;
   
end