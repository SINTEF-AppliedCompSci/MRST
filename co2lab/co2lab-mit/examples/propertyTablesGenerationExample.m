%% Example CO2-brine props for simulation
% Here, we provide a complete example in which we generate our pvt
% properties, put them in a .DATA file, and then create a fluid object for
% simulation with the ad-blackoil module.
% We plot the resulting density and viscosity, as a function of pressure
% and amount of dissolved CO2 in brine.
%
% NOTE: Documentation on the functionality itself is provided within
% pvtBrineWithCO2BlackOil.
%

%% Generate properties
% Properties are based on functional relationships taken from the
% literature. When tables are to be generated, these functions are
% occasionally evaluated outside the strict range of validity suggested in
% the papers to avoid somewhat arbitrary extrapolation when the tables are
% used in a simulation. A warning will be issued when this occurs, but this
% can safely be ignored.
T = {'C', 80};                          % temperature
P = {'MPa','mMn',[0.1, 60, 30]};        % pressure range
S = {'ppm', 'brine', 10^5};             % salinity
saltVar = true;                         % see pvtBrineWithCO2BlackOil
vapH2O = false;                         % "
unsatVals = true;                       % "
figs = false;                           % no figs
[t, rho_co2_s, rho_brine_s] = pvtBrineWithCO2BlackOil(T, P, S, saltVar, ...
                                                      vapH2O, unsatVals, figs);
                                 
%% Load .DATA input file and generate MRST fluid object fior simulation
% Note that we have already copied the output into the .DATA input file that 
% we are going to use, you can check the numbers. We use keywords PVDG
% (dry gas) for the CO2-rich phase, and PVTO (live oil) for the brine which
% can incorporate dissolved CO2.
% Note also that, in the .DATA input file, we remove the last two Rs
% values, since we need to provide data for undersaturated "oil" for the 
% highest Rs (see keyword PVTO in ECLIPSE reference manual, available at 
% http://www.ipt.ntnu.no/~kleppe/TPG4150/EclipseReferenceManual.pdf as of
% July 2022)
% Finally, the brine and CO2 density at surface conditons (defined as
% ECLIPSE standard) must be updated each time, see 2nd and 3rd output in
% pvtBrineWithCO2BlackOil and .DATA input file.

% Load required MRST modules
mrstModule add ad-blackoil deckformat ad-core ad-props
mrstVerbose on
% Load deck from .DATA file
directory = fullfile(getDatasetPath('co2labmit'), 'input_DATA_files');
if unsatVals
    fn = fullfile(directory, 'example_co2brine_3reg_unsatVals.DATA');
else
    fn = fullfile(directory, 'example_co2brine_3reg.DATA');  % less accurate brine properties
end
deck = convertDeckUnits(readEclipseDeck(fn));
% Get fluid object with properties, to be used in simulation
fluid = initDeckADIFluid(deck);     

%% Plot density and viscosity
% Compute values from blackoil-type input properties, and compare with
% values directly reported in the output table from pvt function.
np          = 100;
p_val       = linspace(1,600,np)'*barsa;     % MRST always SI units (Pa)
rho_co2     = fluid.rhoGS*fluid.bG(p_val);   % Dry gas (no water in gas phase)
mu_co2      = fluid.muG(p_val);
rss_val     = fluid.rsSat(p_val);            % Live oil (aqueous phase with dissolved CO2)
rho_b_sat   = fluid.bO(p_val,rss_val,true(np,1)) .* ...
                      (rss_val.*fluid.rhoGS + fluid.rhoOS);
rho_b       = fluid.rhoOS*fluid.bO(p_val,zeros(np,1),false(np,1));
%rho_b       = fluid.rhoOS*interp1(t.P_bar, 1./t.aq.B_b_m3Sm3, p_val/barsa);
mu_b_sat    = fluid.muO(p_val,rss_val,true(np,1));
mu_b        = fluid.muO(p_val,zeros(np,1),false(np,1));

% CO2
latx = {'Interpreter','latex'};
col_co2 = [80, 0, 0; 220, 0, 0]/255;
col_b = [0, 0, 80; 0, 80, 80; 150, 150, 255; 0, 255, 255]/255;
f1 = figure(1);
tiledlayout(1,4,"TileSpacing","compact","Padding","compact")
nexttile(1)
plot(p_val/barsa, rho_co2, '-', 'color', col_co2(1,:), 'linewidth', 3.5); 
hold on
plot(t.P_bar, t.gas.rho_g_kgm3, '-', 'color', col_co2(2,:), 'linewidth', 1.25)
grid on
xlabel('$p$ [bar]', latx{:}, 'fontsize', 12)
ylabel('$\rho_\mathrm{g}$ [kg/m$^3$]', latx{:}, 'fontsize', 12)
xlim([0 P{3}(2)*10])
ylim([0 1200]),
title(['CO$_2$, T=' num2str(T{2}) T{1}], latx{:}, 'fontsize', 14)
legend('fluid object', '\texttt{pvtBrine} fcn', latx{:}, 'fontsize', 12, 'location', 'northwest')
nexttile(2)
plot(p_val/barsa, mu_co2*1e3, '-', 'color', col_co2(1,:), 'linewidth', 3.5); 
hold on
plot(t.P_bar, t.gas.mu_g_cP, '-', 'color', col_co2(2,:), ...
     'linewidth', 1.25)
grid on
xlabel('$p$ [bar]', latx{:}, 'fontsize', 12)
ylabel('$\mu_\mathrm{g}$ [cP]', latx{:}, 'fontsize', 12)
xlim([0 P{3}(2)*10])
ylim([0.01 1])
set(gca,'YScale','log')
%title('CO$_2$ viscosity', latx{:}, 'fontsize', 14)

% Brine 
nexttile(3)
hold on
plot(p_val/barsa, rho_b, '-', 'color', col_b(1,:), 'linewidth', 3.5); 
plot(p_val/barsa, rho_b_sat, '-', 'color', col_b(2,:), 'linewidth', 3.5) 
plot(t.P_bar, t.aq.rho_b_kgm3,'-', 'color', col_b(3,:), 'linewidth', 1.25)
plot(t.P_bar, t.aq.rho_aq_kgm3,'-', 'color', col_b(4,:), 'linewidth', 1.25); 
hold off
grid on
xlabel('$p$ [bar]', latx{:}, 'fontsize', 12)
ylabel('$\rho_\mathrm{b}$ [kg/m$^3$]', latx{:}, 'fontsize', 12)
xlim([0 P{3}(2)*10])
ylim([1000 1100]),
title(['Brine, T=' num2str(T{2}) T{1} ', S=', ...
      num2str(S{3}, '%.0e') S{1}], latx{:}, 'fontsize', 14)
legend('no CO$_2$', 'sat. CO$_2$', 'no CO$_2$, fcn', ...
       'sat. CO$_2$, fcn', latx{:}, 'fontsize', 12, 'location', 'southeast')
nexttile(4)
hold on
plot(p_val/barsa, mu_b*1e3, '-', 'color', col_b(1,:), 'linewidth', 3.5); 
plot(p_val/barsa, mu_b_sat*1e3, '-', 'color', col_b(2,:), 'linewidth', 3.5)
plot(t.P_bar, t.aq.mu_b_cP, '-', 'color', col_b(3,:), 'linewidth', 1.25); 
plot(t.P_bar, t.aq.mu_aq_cP, '-', 'color', col_b(4,:), 'linewidth', 1.25);
hold off
grid on
xlabel('$p$ [bar]', latx{:}, 'fontsize', 12)
ylabel('$\mu_\mathrm{b}$ [cP]', latx{:}, 'fontsize', 12)
xlim([0 P{3}(2)*10])
ylim([0.3 0.6])
%set(gca,'YScale','log')
%title('Brine viscosity', latx{:}, 'fontsize', 14)
set(f1, 'Position', [200, 200, 800, 250])