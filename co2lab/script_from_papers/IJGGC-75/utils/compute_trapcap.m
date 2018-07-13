function [H1, strap, btrap_res, btrap_dis, p] = compute_trapcap(Gt, ta, rock, seainfo, surface_pressure, varargin)
% local helper function adapted from exploreCapacity

  opt.press_deviation = 0;  % deviation (pos or neg percentage) which gets added to hydrostatic condition
  opt.caprockTemp = [];     % if empty, will be computed (in Kelvin).
  opt = merge_options(opt, varargin{:});
  
  poro = rock.poro;
  ntg = ones(Gt.cells.num,1);
  if isfield(rock,'ntg')
     ntg = rock.ntg; 
  end

  % Computing co2 density as function of reservoir conditions
  p = computeHydrostaticPressure(Gt, seainfo.water_density, surface_pressure); % Pa
  pv = poro .* Gt.cells.volumes .* Gt.cells.H .* ntg; % m^3
  p = p + sum(p.*pv).*opt.press_deviation/(100.*sum(pv));
  if isempty(opt.caprockTemp)
      t = computeCaprockTemperature( Gt, seainfo.seafloor_temp, ...
                                 seainfo.seafloor_depth, ...
                                 seainfo.temp_gradient ); % Celsius
      t = t + 273.15; % Kelvin
  else
      t = opt.caprockTemp; 
  end
  co2 = CO2props('sharp_phase_boundary', true);
  co2_rho = co2.rho(p,t);

  % Computing structural trap heights (H1) for each cell
  trapcells     = find(ta.traps);
  H1            = zeros(Gt.cells.num, 1);
  if ~isempty(trapcells)
     H1(trapcells) = ta.trap_z(ta.traps(trapcells)) - Gt.cells.z(trapcells);
  end
  H1=min(H1,Gt.cells.H);
  H2 = Gt.cells.H - H1;
  assert(all(H1<=Gt.cells.H));

  % Computing total trapping volume in structural traps (dissolved and
  % structurally trapped
  strap_pvol_tot       = Gt.cells.volumes .* H1 .* poro .* ntg;
  strap_pvol_co2_plume = strap_pvol_tot .* (1 - seainfo.res_sat_wat);
  strap_pvol_co2_diss  = strap_pvol_tot .* seainfo.res_sat_wat .* seainfo.dis_max;

  strap = strap_pvol_co2_plume .* co2_rho + ...
          strap_pvol_co2_diss  .* seainfo.rhoCref;

  % Computing total trapping volume below structural traps (dissolved and
  % residually trapped
  btrap_pvol_tot          = Gt.cells.volumes .* H2 .* poro;
  btrap_pvol_co2_residual = btrap_pvol_tot .* seainfo.res_sat_co2;
  btrap_pvol_co2_dissol   = btrap_pvol_tot .* (1-seainfo.res_sat_co2) .* seainfo.dis_max;

  btrap_res = btrap_pvol_co2_residual .* co2_rho;
  btrap_dis = btrap_pvol_co2_dissol   .* seainfo.rhoCref;
end