function [ hfig, hax ] = plotProfilesAtGivenTime(givenTime, Gt, states, initState, fluid, model, caprock_temperature)

% givenTime can be either 'final' (type char) or a specified time (type
% double) in seconds after simulation start

if strcmpi(givenTime,'final')
    % meaningful profiles
    press       = states{end}.pressure;
    pressDiffFromHydrostatic = press - initState.pressure;
    densityCO2  = fluid.rhoG(states{end}.pressure);  % fluid.rhoG is function handle to get CO2 density
    satCO2      = states{end}.s(:,2);
    massCO2     = model.rock.poro.*model.G.cells.volumes.* model.G.cells.H.*satCO2.*densityCO2; % kg
    
elseif isa(givenTime,'double')
    [rti,~] = find(sim_report.ReservoirTime==givenTime);
    
    press       = states{rti}.pressure;
    pressDiffFromHydrostatic = press - initState.pressure;
    densityCO2  = fluid.rhoG(states{rti}.pressure);  % fluid.rhoG is function handle to get CO2 density
    satCO2      = states{rti}.s(:,2);
    massCO2     = model.rock.poro.*model.G.cells.volumes.* model.G.cells.H.*satCO2.*densityCO2; % kg
end


bf = boundaryFaces(Gt);

figure; set(gcf,'Position',[1 1 3000 400])
hold on

subplot(1,6,1)
plotCellData(Gt, caprock_temperature, 'EdgeColor','none')
title('Caprock Temperature'); axis off equal tight
[ ~ ] = setColorbarHandle( gcf, 'LabelName', 'Kelvin', 'fontSize', 18 );


subplot(1,6,2)
plotCellData(Gt, press, 'EdgeColor','none')
title('Pressure'); axis off equal tight
[ ~ ] = setColorbarHandle( gcf, 'LabelName', 'Pascals', 'fontSize', 18 );


subplot(1,6,3)
plotCellData(Gt, pressDiffFromHydrostatic, 'EdgeColor','none')
title({'Pressure diff. from hydrostatic';'i.e., the initial condition'}); axis off equal tight
[ ~ ] = setColorbarHandle( gcf, 'LabelName', 'Pascals', 'fontSize', 18 );


subplot(1,6,4)
plotCellData(Gt, densityCO2, 'EdgeColor','none')
title('CO2 density at Caprock'); axis off equal tight
[ ~ ] = setColorbarHandle( gcf, 'LabelName', 'kg/m^3', 'fontSize', 18 );


subplot(1,6,5)
%plotGrid(Gt, 'FaceColor', 'white', 'EdgeAlpha', 0)
plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
plotCellData(Gt, satCO2, satCO2>(0.1/100), 'EdgeColor','none') %satCO2~=0)
title('Saturation of CO2'); axis off equal tight
[ ~ ] = setColorbarHandle( gcf, 'LabelName', 'Pore Volume CO2 Saturation', 'fontSize', 18 );


subplot(1,6,6)
%plotGrid(Gt, 'FaceColor', 'white', 'EdgeAlpha', 0)
plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
plotCellData(Gt, massCO2/1e9, satCO2>(0.1/100), 'EdgeColor','none') % only plot plume that has sat > 0.1 percent
title('Mass of CO2'); axis off equal tight
[ ~ ] = setColorbarHandle( gcf, 'LabelName', 'Mt', 'fontSize', 18 );


hfig = gcf;
hax  = gca;

end

