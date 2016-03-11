function compare_p_impact_on_retaining_cap( Gt, optim, other )
%Study impact of pressure on max retaining capacity breakdown and trapping achieved.

% DESCRIPTION
%   Pressure impacts CO2 density, and thus impacts the amount of CO2 that
%   may be structurally trapped. The max retaining capacity is computed for
%   hydrostatic pressure (initState) and then for the final pressure field
%   of an injection simulation. The amount of trapping achieved is then
%   computed relative to each pressure scenario.

    % co2 density as a function of hydrostatic pressure:
    %   other.fluid.rhoG(other.initState.pressure)

    % co2 density as a function of final pressure state:
    %   other.fluid.rhoG(optim.states{end}.pressure)
    % or
    %   other.fluid.rhoGS .* other.fluid.bG(optim.states{end}.pressure)


    if ~isfield(other,'fluid')
        other.fluid = reconstructFluidStructure( Gt, other.rock, other.opt, other.dh );
    end
    
    seainfo = getSeaInfo(other.opt.modelname, other.opt.refRhoCO2);
    
    % deviation of final pressure from hydrostatic
    press_dev = (other.fluid.rhoG(optim.states{end}.pressure) - ...
            other.fluid.rhoG(other.initState.pressure)) ...
            ./other.fluid.rhoG(other.initState.pressure) * 100;

    % Max retaining capacity given hydrostatic pressure:
    % (could pass in initState.pressure as a pressure field, or let it be
    % computed internally)
    fmCap = getTrappingInfo(Gt, other.rock, seainfo, ...
              'mapPlotOn',false, 'surf_press',other.opt.surface_pressure);
    
    % Max retaining capacity given final pressure state:
    % (could pass in final pressure state as a pressure field, or set
    % seainfo.press_dev = press_dev as computed above)
    %seainfo.press_deviation = press_dev;
    fmCap = getTrappingInfo(Gt, other.rock, seainfo, ...
              'mapPlotOn',false, 'press_field',optim.states{end}.pressure);
    
end

