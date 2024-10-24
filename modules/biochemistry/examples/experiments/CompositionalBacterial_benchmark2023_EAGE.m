%% MRST Simulation for Hydrogen Storage with Bacterial Growth Model
% Description: This script uses MRST to model gas injection into a 3D porous medium,
% incorporating compositional fluid properties, bacterial mono modal.
% We consider a liquid phase (W) and a gas (G) phase, 4 components 
% ('H2O','H2','CO2','CH4') and The microbial activity of 
%a archaea.
%This test case comes from a Benchmark in EAGE 2023
% Clear workspace and initialize MRST modules
clear; clc;
mrstModule add biochemistry compositional ad-blackoil ad-core ad-props mrst-gui
gravity reset on 
biochemistrymodel=true; 
writedatafile=false;
useHandler=false;%true;% 
compare_bact=true;%false;%
%% ============Grid and Rock Properties=====================
% Define grid dimensions and physical dimensions
%[nx, ny, nz] = deal(61,61,10);  % Grid cells in x, y, z directions
[nx, ny, nz] = deal(31,31,10);  % Grid cells in x, y, z directions
[Lx,Ly,Lz] = deal(1525,1525,50);         % Physical dimensions in meters
dims = [nx, ny, nz];
pdims = [Lx, Ly, Lz];


% Create grid and shift vertically by reservoir depth
G = cartGrid(dims, pdims);
depth_res = 1026;   % Reservoir depth in meters
G.nodes.coords(:, 3) = G.nodes.coords(:, 3) + depth_res;
G = computeGeometry(G);

% Define rock properties
K=[100, 100, 10].*milli*darcy;
rock = makeRock(G, K, 0.2);  % Default permeability and porosity


%% Fluid Properties Initialization
% Define compositional fluid model (with CoolProp library support)
compFluid = TableCompositionalMixture({'Water', 'Hydrogen', 'CarbonDioxide', 'Methane'}, ...
                                      {'H2O', 'H2', 'CO2', 'C1'});

% Fluid density and viscosity (kg/m^3 and cP)
[rhow, rhog] = deal(999.7 * kilogram / meter^3, 0.0817 * kilogram / meter^3);
[viscow, viscog] = deal(1.3059 * centi * poise, 0.01763 * centi * poise);

% Compressibility (per bar)
[cfw, cfg] = deal(4.1015e-5/ barsa, 8.0e-3 / barsa);

% Relative permeability and initial saturations
%[srw, src] = deal(0.0, 0.0);
[srw, src] = deal(0.2, 0.05);
P0=100 * barsa;%106 * barsa;
fluid = initSimpleADIFluid('phases', 'OG', 'mu', [viscow, viscog], ...
                           'rho', [rhow, rhog], 'pRef', P0, ...
                           'c', [cfw, cfg], 'n', [2, 2], 'smin', [srw, src]);

% Capillary pressure function
Pe = 0.1 * barsa;
pcOG = @(sw) Pe * sw.^(-1/2);
fluid.pcOG = @(sg) pcOG(max((1 - sg - srw) / (1 - srw), 1e-5));

%% Simulation Parameters
% Set total time, pore volume, and injection rate
niter=140;
TotalTime = niter*day;
rate = 1e6*meter^3/day; 


%% Time Stepping and Schedule
% Define schedule and solver
nls = NonLinearSolver('useRelaxation', true);
deltaT =rampupTimesteps(TotalTime, 1*day, 0);
schedule = simpleSchedule(deltaT);
nj1=90;nj2=100;nj3=130;
schedule.step.control(1:nj1)=1;
schedule.step.control(nj1+1:nj2)=2;
schedule.step.control(nj2+1:nj3)=3;
schedule.step.control(nj3+1:end)=4;

%% Wells and Boundary Conditions
% Initialize wells
W1 = [];
W2 = [];
W3 = [];
W4 = [];
tmp = cell(4,1);
n1=floor(0.5*nx)+1; n2=floor(0.5*nx)+1;
schedule.control = struct('W',tmp);
cellInd=zeros(nz-1,1);
for k=2:nz
    cellInd(k-1)=(k-1)*nx*ny+(n2-1)*nx+n1;
end
% Injection well parameters
W1 = verticalWell(W1, G, rock, n1, n2, 1:nz-1, 'comp_i', [0, 1], 'Radius', 0.5, ...
                 'name', 'Injector', 'type', 'rate', 'Val', rate, 'sign', 1);
W1(1).components = [0.0, 0.95,  0.05, 0.0];  % H2-rich injection   {'H2O', 'H2', 'CO2', 'C1'});

%Idle period
W2 = verticalWell(W2, G, rock, n1, n2, 1:nz-1, 'compi', [0, 1], 'Radius', 0.5, ...
                 'name', 'Rest', 'type', 'rate', 'Val', 0.0, 'sign', 1);
W2(1).components = [0.0, 0.95,  0.05, 0.0];  % rest period

%production
%Pwell=66*barsa;
%W3 = verticalWell(W3, G, rock, n1, n2, 1:nz-1, 'comp_i', [0, 1], 'Radius', 0.5, ...
 %                 'name', 'Prod', 'type', 'bhp', 'Val', Pwell, 'sign', -1);
W3 = verticalWell(W3, G, rock, n1, n2, 1:nz-1, 'compi', [0, 1], 'Radius', 0.5, ...
                  'name', 'Prod', 'type', 'rate', 'Val', -rate, 'sign', -1);
W3(1).components = [0.0, 0.95,  0.05, 0.0];  %production

%Idle period
W4 = verticalWell(W4, G, rock, n1, n2, 1:nz-1, 'compi', [0, 1], 'Radius', 0.5, ...
                 'name', 'Rest', 'type', 'rate', 'Val', 0.0, 'sign', 1);
W4(1).components = [0.0, 0.95,  0.05, 0.0];  % rest period


schedule.control(1).W = W1;
schedule.control(2).W = W2;
schedule.control(3).W = W3;
schedule.control(4).W = W4;

%% Model Setup: Compositional Model with Bacterial Growth
if biochemistrymodel
    eosname='sw';
    eosmodel =SoreideWhitsonEquationOfStateModel(G, compFluid,eosname);
    diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);
    mex_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'useMex', true, 'rowMajor', true);

    arg = {G, rock, fluid, compFluid,true,diagonal_backend,...
        'water', false, 'oil', true, 'gas', true,'bacteriamodel', true,...
        'bDiffusionEffect', false,'moleculardiffusion',false,...
        'liquidPhase', 'O', 'vaporPhase', 'G'};
    model = BiochemistryModel(arg{:});
    model.outputFluxes = false;
    % model.Psigrowthmax=5.e-5;
    model.b_bact= 1.e-10;
    model.Y_H2 = 1.010e12;
    model.alphaH2 = 1.1e-6;
    model.alphaCO2 = 3.2e-5;
    model.Psigrowthmax = 1.7e-4;
    model.nbactMax = 6.88e11; 
    model.EOSModel.msalt=0;
else
    eosname='pr';
    arg = {G, rock, fluid, compFluid, 'water', false, 'oil', true, 'gas', true, ...
      'liquidPhase', 'O', 'vaporPhase', 'G'};
    model = GenericOverallCompositionModel(arg{:});
end



%% Initial Conditions
% Temperature and initial saturations
T0 = 40+273.15;                % Initial temperature (K)
s0 = [0.2, 0.8];           % Initial saturations (Sw,Sg)
%z0 = [0.2, 0.0, 0.0, 0.8];  % Initial composition: H2O, H2, CO2, CH4
z0 = [0.82, 1.0e-5, 1.0e-5, 0.18-2.0e-5];  % Initial composition: H2O, H2, CO2, CH4
Phydro0=rhow*norm(gravity).*G.cells.centroids(:,3);
% Initialize state with bacterial concentration

if biochemistrymodel
    if model.bacteriamodel
        nbact0 = 1.0e11;%1e6;  
        state0 = initCompositionalStateBacteria(model, P0, T0, s0, ...
            z0, nbact0,eosmodel);
    else
        state0 = initCompositionalState(model, Phydro0, T0, s0, z0);
    end

else
    state0 = initCompositionalState(model, Phydro0, T0, s0, z0);
end


%% Run simulation
mrstModule add mpfa
model_mpfa = setMPFADiscretization(model);
%[wellSols,states,report]= simulateScheduleAD(state0, model, schedule, 'nonlinearsolver', nls);
if useHandler 
    dir='/home/sdelage2/PROJETS/gdr_h2/MRST2024/MRST/output';
    diroutput='Benchmark2023AEGE_NOBACT';
    handler = ResultHandler('writeToDisk', true,'dataDirectory',dir,...
        'dataFolder', diroutput);
    [wellSols,states,report]= simulateScheduleAD(state0, model, schedule,...
        'nonlinearsolver',nls,'outputHandler', handler);
else
    [wellSols,states,report]= simulateScheduleAD(state0, model_mpfa, schedule, 'nonlinearsolver', nls);
end  


%% Plotting Results and estimate bacteria effects
namecp = model.EOSModel.getComponentNames();
indH2=find(strcmp(namecp,'H2'));
indCO2= find(strcmp(namecp,'CO2'));
indCH4= find(strcmp(namecp,'C1'));
nT=numel(states);
xH2=zeros(nT,1);
yH2=zeros(nT,1);
xCO2= zeros(nT,1);
yCO2= zeros(nT,1);
yCH4= zeros(nT,1);
pressure=zeros(nT,1);

totMassH2= zeros(nT,1);
totMassCO2= zeros(nT,1);
totMassCH4= zeros(nT,1);
FractionMassCO2= zeros(nT,1);
FractionMassH2= zeros(nT,1);
FractionMassCH4= zeros(nT,1);
totMassComp= zeros(nT,1);
ncomp=model.EOSModel.getNumberOfComponents;

for i = 1:nT
    xH2(i)=max(states{i}.x(:,indH2));
    yH2(i)=max(states{i}.y(:,indH2));
    yCH4(i)=max(states{i}.y(:,indCH4));
    xCO2(i)=max(states{i}.x(:,indCO2));
    yCO2(i)=max(states{i}.y(:,indCO2));
    pressure(i)=mean(states{i}.pressure(:));
end

for i = 1:nT
    for j=1:ncomp
        totMassComp(i)=totMassComp(i)+sum(states{i}.FlowProps.ComponentTotalMass{j});
    end
    totMassH2(i)=sum(states{i}.FlowProps.ComponentTotalMass{indH2});
    totMassCO2(i)=sum(states{i}.FlowProps.ComponentTotalMass{indCO2});
    totMassCH4(i)=sum(states{i}.FlowProps.ComponentTotalMass{indCH4});
    FractionMassH2(i)=totMassH2(i)/totMassComp(i);
    FractionMassCO2(i)=totMassCO2(i)/totMassComp(i);
    FractionMassCH4(i)=totMassCH4(i)/totMassComp(i);
end
%% Calculate H2 production efficiency
mH2_injected=totMassH2(nj1)-totMassH2(1);
mH2_produced=totMassH2(nj2+1)-totMassH2(nj3);
Efficiency_H2=(mH2_produced/mH2_injected) * 100;
fprintf('H2 Production Efficiency : %.2f%%\n', Efficiency_H2);

%% Compare case without bacteria and with bacteria
if compare_bact
    %Extraction data states_nobact de handler
    dir='/home/sdelage2/PROJETS/gdr_h2/MRST2024/MRST/output';
    handler1 = ResultHandler('dataDirectory',dir,'dataFolder','Benchmark2023AEGE_NOBACT');
    m = handler1.numelData();
    states_nobact = cell(m, 1);
    for i = 1:m
        states_nobact{i} = handler1{i};
    end
    totMassH2_nobact= zeros(nT,1);
    totMassCO2_nobact= zeros(nT,1);
    totMassCH4_nobact= zeros(nT,1);
    totMassComp_nobact= zeros(nT,1);
    FractionMassCO2_nobact= zeros(nT,1);
    FractionMassH2_nobact= zeros(nT,1);
    FractionMassCH4_nobact= zeros(nT,1);
    pressure_nobact=zeros(nT,1);
    xH2_nobact=zeros(nT,1);
    xCO2_nobact= zeros(nT,1);
    yH2_nobact=zeros(nT,1);
    yCO2_nobact= zeros(nT,1);
    yCH4_nobact= zeros(nT,1);

    %total mass
    for i = 1:nT
        for j=1:ncomp
        totMassComp_nobact(i)=totMassComp_nobact(i)+sum(states_nobact{i}.FlowProps.ComponentTotalMass{j});
        end
        totMassH2_nobact(i)=sum(states_nobact{i}.FlowProps.ComponentTotalMass{indH2});
        totMassCO2_nobact(i)=sum(states_nobact{i}.FlowProps.ComponentTotalMass{indCO2});
        totMassCH4_nobact(i)=sum(states_nobact{i}.FlowProps.ComponentTotalMass{indCH4});
        FractionMassH2_nobact(i)=totMassH2_nobact(i)/totMassComp(i);
        FractionMassCO2_nobact(i)=totMassCO2_nobact(i)/totMassComp(i);
        FractionMassCH4_nobact(i)=totMassCH4_nobact(i)/totMassComp(i);
        pressure_nobact(i)=mean(states_nobact{i}.pressure(:));
        xH2_nobact(i)=max(states_nobact{i}.x(:,indH2));
        xCO2_nobact(i)=max(states_nobact{i}.x(:,indCO2));
        yH2_nobact(i)=max(states_nobact{i}.y(:,indH2));
        yCO2_nobact(i)=max(states_nobact{i}.y(:,indCO2));
        yCH4_nobact(i)=max(states_nobact{i}.y(:,indCH4));
    end

    %% Calculate H2 production efficiency 
    mH2_injected_nobact=totMassH2_nobact(nj1)-totMassH2_nobact(1);
    mH2_produced_nobact=totMassH2_nobact(nj2+1)-totMassH2_nobact(nj3);
    Efficiency_H2_nobact=(mH2_produced_nobact/mH2_injected_nobact) * 100;
    %% Display H2 production efficiency
    fprintf('H2 Production Efficiency without bacteria: %.2f%%\n', Efficiency_H2_nobact);
    fprintf('H2 Production Efficiency with bacteria: %.2f%%\n', Efficiency_H2);
  
    %% Calculate percentage of H2 loss
    H2_loss_percentage = (abs(totMassH2_nobact-totMassH2)./totMassH2_nobact) * 100;
    %% Calculate percentage of CO2 loss
    CO2_loss_percentage = (abs(totMassCO2_nobact-totMassCO2)./totMassCO2_nobact) * 100;
    %% Calculate percentage of CH4 production
    CH4_loss_percentage = (abs(totMassCH4_nobact-totMassCH4)./totMassCH4_nobact) * 100;

    %% Display final H2 loss
    fprintf('Total H2 loss due to bacterial effects: %.2f%%\n', H2_loss_percentage(end));
    fprintf('Total CO2 loss due to bacterial effects: %.2f%%\n', CO2_loss_percentage(end));
    fprintf('Total CH4 production due to bacterial effects: %.2f%%\n', CH4_loss_percentage(end));

    %% plot mass
    figure(1); clf; 
    plot(1:nT,FractionMassH2_nobact,'b')
    hold on
    plot(1:nT,FractionMassH2,'k-')
    title('H2 Mass fractions')
    xlabel('Time (days)')
    ylabel('mass fraction')
    legend('YH2 nobact','YH2 with bact')

    figure(2); clf; 
    plot(1:nT,totMassH2_nobact,'b')
    hold on
    plot(1:nT,totMassCO2_nobact,'k-')
    title('Mass ')
    xlabel('Time (days)')
    ylabel('mass fraction')
    legend('mH2','mCO2')

   

f1=figure('Name','Pressure_compar_bact','NumberTitle','off');
f1.Position(3:4) = [900 700];
plot(1:nT,pressure./1e5,'b','MarkerSize',7,'LineWidth',2)
hold on
plot(1:nT,pressure_nobact./1e5,'r--','MarkerSize',8,'LineWidth',2)

title('Mean pressure in the reservoir','FontSize',14,'FontWeight','bold','Color','k')
xlabel({'time','(days)'},'FontWeight','bold','Color','k')
ylabel({'pressure','(bar)'},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 14; 

legend({'pressure, with archae','pressure, no archae'},...
    'FontSize',13,'TextColor','black',...
    'Location','best')
%xlim([1 nT])
%ylim([min(min_xliqH2)-1e-4 max(max_xliqH2)+1e-4])


f2=figure('Name','maxyH2_compar_bact','NumberTitle','off');
f2.Position(3:4) = [900 700];
plot(0:nT,[0;yH2],'b','MarkerSize',7,'LineWidth',2)
hold on
plot(0:nT,[0;yH2_nobact],'r--','MarkerSize',8,'LineWidth',2)

title('Maximum H_2 Molar fractions in gaz','FontSize',14,'FontWeight','bold','Color','k')
xlabel({'time','(days)'},'FontWeight','bold','Color','k')
ylabel({'H_2 molar fraction',''},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 14; 

legend({'X_{g,H_2}, with archae','X_{g,H_2}, no archae'},...
    'FontSize',13,'TextColor','black',...
    'Location','best')


f20=figure('Name','maxxH2_compar_bact','NumberTitle','off');
f20.Position(3:4) = [900 700];
plot(0:nT,[0;xH2],'b','MarkerSize',7,'LineWidth',2)
hold on
plot(0:nT,[0;xH2_nobact],'r--','MarkerSize',8,'LineWidth',2)

title('Maximum H_2 Molar fractions in liquid','FontSize',14,'FontWeight','bold','Color','k')
xlabel({'time','(days)'},'FontWeight','bold','Color','k')
ylabel({'H_2 molar fraction',''},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 14; 

legend({'X_{a,H_2}, with archae','X_{a,H_2}, no archae'},...
    'FontSize',13,'TextColor','black',...
    'Location','best')

f21=figure('Name','maxxCO2_compar_bact','NumberTitle','off');
f21.Position(3:4) = [900 700];
plot(0:nT,[0;xCO2],'b','MarkerSize',7,'LineWidth',2)
hold on
plot(0:nT,[0;xCO2_nobact],'r--','MarkerSize',8,'LineWidth',2)

title('Maximum CO_2 Molar fractions in liquid','FontSize',14,'FontWeight','bold','Color','k')
xlabel({'time','(days)'},'FontWeight','bold','Color','k')
ylabel({'CO_2 molar fraction',''},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 14; 

legend({'X_{a,CO_2}, with archae','X_{a,CO_2}, no archae'},...
    'FontSize',13,'TextColor','black',...
    'Location','best')

f3=figure('Name','maxyCO2_compar_bact','NumberTitle','off');
f3.Position(3:4) = [900 700];
plot(0:nT,[0;yCO2],'b','MarkerSize',7,'LineWidth',2)
hold on
plot(0:nT,[0;yCO2_nobact],'r--','MarkerSize',8,'LineWidth',2)

title('Maximum CO_2 Molar fractions in gaz','FontSize',14,'FontWeight','bold','Color','k')
xlabel({'time','(days)'},'FontWeight','bold','Color','k')
ylabel({'CO_2 molar fraction',''},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 14; 

legend({'X_{g,CO_2}, with archae','X_{g,CO_2}, no archae'},...
    'FontSize',13,'TextColor','black',...
    'Location','best')


f4=figure('Name','fractionmass_compar_bact','NumberTitle','off');
f4.Position(3:4) = [900 700];
plot(0:nT,[0;FractionMassH2],'b','MarkerSize',7,'LineWidth',2)
hold on
plot(0:nT,[0;FractionMassH2_nobact],'r--','MarkerSize',8,'LineWidth',2)

title('Total H_2 mass fraction ','FontSize',14,'FontWeight','bold','Color','k')
xlabel({'time','(days)'},'FontWeight','bold','Color','k')
ylabel({'H_2 mass fraction',''},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 14; 

legend({'X_{H_2}, with archae','X_{H_2}, no archae'},...
    'FontSize',13,'TextColor','black',...
    'Location','best')



end %end comparebact


figure(3); clf; 
plot(1:nT,totMassH2,'b')
hold on
plot(1:nT,totMassCO2,'k-')
title('Mass ')
xlabel('Time (days)')
ylabel('mass of H_2 and CO_2')
legend('mH2','mCO2')

figure(4); clf; 
plot(1:nT,FractionMassH2,'b')
hold on
plot(1:nT,FractionMassCO2,'k-')
title('Mass fractions with bacteria')
xlabel('Time (days)')
ylabel('mass fraction')
legend('YH2','YCO2')



figure(5); clf; 
plot(1:nT,pressure./1e5,'b')
title('pressure')
xlabel('Time (days)')
ylabel('pressure (bar)')
legend('pressure')

figure(6); clf; 
plot(1:nT,yH2,'b')
hold on
plot(1:nT,yCO2,'k-')
title('Maximum Molar fractions in gas')
xlabel('Time (days)')
ylabel('molar fraction')
legend('yH2','yCO2')

figure(7); clf; 
plot(1:nT,xH2,'b')
hold on
plot(1:nT,xCO2,'k-')
title('Maximum Molar fractions in Liquid')
xlabel('Time (days)')
ylabel('molar fraction')
legend('xH2','xCO2')


figure(7); clf; 
plot(0:nT,[0;xH2/max(xH2)],'b')
hold on
plot(0:nT,[0;xCO2/max(xCO2)],'k-')
title('Maximum Molar fractions in Liquid')
xlabel('Time (days)')
ylabel('molar fraction')
legend('xH2/xH2_{max}','xCO2/xCO2_{max}')

figure(7); clf; 
plot(0:nT,[0;yH2/max(yH2)],'b')
hold on
plot(0:nT,[0;yCO2/max(yCO2)],'k-')
title('Maximum Molar fractions in gaz')
xlabel('Time (days)')
ylabel('molar fraction')
legend('yH2/yH2_{max}','yCO2/yCO2_{max}')





if biochemistrymodel && model.bacteriamodel
    nbacteria= zeros(nT,1);
    pv=model.operators.pv;
    ncells=G.cells.num;
    for i = 1:nT
       Swi = states{i}.s(:,1);
       Swpvi=Swi.*pv;
       Swpv=sum(Swpvi);
       nbacteria(i)=sum(states{i}.nbact.*Swpvi)/Swpv;
       %nbacteria(i)=max(states{i}.nbact);
    end

f31=figure('Name','nbacteria','NumberTitle','off');
f31.Position(3:4) = [900 700];
plot(0:nT,[nbact0;nbacteria],'b','MarkerSize',7,'LineWidth',2)

title('Methanogenic Archae','FontSize',14,'FontWeight','bold','Color','k')
xlabel({'time','(days)'},'FontWeight','bold','Color','k')
ylabel({'n_{archae}','(m^{-3})'},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 14; 

legend({'n_{archae}'},...
    'FontSize',13,'TextColor','black',...
    'Location','best')



f11=figure('Name','nbacteria_t140','NumberTitle','off');
f11.Position(3:4) = [900 700];
plotCellData(G,states{140}.nbact);
colorbar; 
axis equal
axis ([0 Lx  0 Ly depth_res depth_res+Lz])
view(0,-90)
title('Methanogenic Archae density, 140 days','FontSize',14,'FontWeight','bold','Color','k')
xlabel({'x','(m)'},'FontWeight','bold','Color','k')
ylabel({'y','(m)'},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 14; 
%legend({'n_{archae}'},...
 %   'FontSize',13,'TextColor','black',...
%    'Location','east')
clim([0 nbact0])


    figure(8); clf; 
    plot(1:nT,nbacteria,'b')
    title('Bacteria density in the liquid phase')
    xlabel('Time (days)')
    ylabel('bacteria density')
    legend('nbacteria')

    figure(9),clf
    for i=1:nT
        clf;
        plotCellData(G,states{i}.nbact);
        colorbar; 
        axis equal
        axis ([0 Lx  0 Ly depth_res depth_res+Lz])
        view(0,-90)
        pause(0.1)   
    end
    title('Archea density')
    legend('nbacteria')

end

if writedatafile
    outputFileName = sprintf('H2_CO2_solubility_%s_withnbact.dat',eosname);
    fileID = fopen(outputFileName, 'wt');

    fprintf(fileID, 'iter H2_molar_fraction  CO2_molar_fraction\n');
    for i = 1:nT
       fprintf(fileID, '%d  %12.8f  %12.8f %12.8f  %12.8f\n', i, xH2(i),xCO2(i), yH2(i),yCO2(i));
    end
    fclose(fileID);

    disp(['Results written to file: ', outputFileName]);
end
