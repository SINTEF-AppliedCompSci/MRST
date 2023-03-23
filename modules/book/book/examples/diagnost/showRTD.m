%% Illustrate the use of residence-time distributions
% In this example, we compute the residence-time distribution for two
% different layers of the SPE 10 benchmark test. We compare distributions
% obtained by tracing a delta pulse through the reservoir and obtained
% directly from the averaged TOF values. 
mrstModule add diagnostics spe10 incomp

%% Base model
% We set up a grid. Later, we will assign petrophysical properties from one
% single layer at the time
cartDims = [  60,  220,  1];
physDims = [1200, 2200, 2*cartDims(end)] .* ft();
G  = cartGrid(cartDims, physDims);
G  = computeGeometry(G);

% Set fluid model: Here, we use a simple fluid model with properties that
% are typical for water. Replace this by a multiphase fluid object if you
% also want to include fluid effects in the calculation
fluid = initSingleFluid('mu', 1*centi*poise, ...
                        'rho', 1014*kilogram/meter^3);

%% Compute RTD for a single well
[T,layer] = deal(10*year, [23, 75]);
for n=1:2
    % --- Set petrophysical data for this particular layer
    % To avoid problems with very small porosity values, we explicitly
    % impose a lower threshold of 1e-4
    rock = getSPE10rock(1:cartDims(1),1:cartDims(2),layer(n));
    rock.poro = max(rock.poro, 1e-4);
    
    % --- Set up well model
    pv  = poreVolume(G, rock);
    W = addWell([], G, rock, 1:60, ...
            'Type', 'rate', 'Val', 3*sum(pv)/T, ...
            'Radius', 0.125, 'Name', 'I', ...
            'InnerProduct', 'ip_tpf', 'Comp_i', 1);
    W = addWell(W, G, rock, (60*219+1):(60*220), ...
            'Type', 'bhp', 'Val', 200*barsa, ...
            'Radius', 0.125, 'Name', 'P', ...
            'InnerProduct', 'ip_tpf', 'Comp_i', 1);
        
    % --- Initiate and solve flow problem 
    rS = initState(G, W, 0, 0.0);
    hT = computeTrans(G, rock);
    rS = incompTPFA(rS, G, hT, fluid, 'wells', W);
    rS.wellSol(1).sign = 1;
    rS.wellSol(2).sign =-1;
    
    
    % --- Compute RTD
    D   = computeTOFandTracer(rS, G, rock, 'wells', W, ...
        'maxTOF', inf,'computeWellTOFs', true, 'firstArrival', true);
    WP  = computeWellPairs(rS, G, rock, W, D);
    rtd = computeRTD(rS, G, pv, D, WP, W, 'nsteps',5000);
    RTD = estimateRTD(pv, D, WP);
    
    % --- Plot RTD
    figure, hold on
    plot(rtd.t/year, rtd.values, '-',  'LineWidth',2);
    plot(RTD.t/year, RTD.values, '--', 'LineWidth',2);
    tfa = min(D.ifa(W(2).cells))/year;      % first arrival
    tm  = RTD.volumes/RTD.allocations/year; % mean time
    vm  = max(RTD.values);
    plot([tfa tfa], [0 vm], ':k', [tm tm], [0 vm], '-k');
    legend('Simulated pulse','From averaged TOF');    
    set(gca,'XLim',[0 10],'YLim',[0 5e-8],'FontSize',14); 
    xlabel('Time [years]')
    hold off
    
    axes('Position',[.6 .4 .3 .35]);
    plotCellData(G,log10(rock.perm(:,1)),'EdgeColor','none');
    view(2); axis equal off, colormap(parula)
    
    % --- F-Phi diagram
    [F, Phi] = computeFandPhi(pv, D.tof);
    [f, phi] = computeFandPhi(rtd, 'sum', true);

    figure
    plot(phi,f, '-', Phi,F, '--', 'LineWidth',2);
    legend(['Simulated pulse: L=',num2str(computeLorenz(f,phi))],...
        ['From averaged TOF: L=', num2str(computeLorenz(F,Phi))], ...
        'Location','SouthEast');
    axis([0 1 0 1]); set(gca,'FontSize',14);
end

%% RTD for two well pairs
% Set parameters describing the wells
wtype    = {'bhp', 'bhp', 'bhp'};
wtarget  = [200,   200,   600] .* barsa();
wrad     = [0.125, 0.125, 0.125] .* meter;
wloc     = [ 10,   50,    25;
            220,   220,    1];
wname    = {'P1', 'P2', 'I'};
sgn      = [-1 -1 1];

[T,layer] = deal(10*year, [23, 75]);
for n=1:2
    
    % --- Set petrophysical data for this particular layer
    % To avoid problems with very small porosity values, we explicitly
    % impose a lower threshold of 1e-4
    rock = getSPE10rock(1:cartDims(1),1:cartDims(2),layer(n));
    rock.poro = max(rock.poro, 1e-4);
    
    % --- Set up well model
    % To ensure that we get the correct well index when updating the
    % petrophysical data, we simply regenerate the well objects.
    W = [];
    for w = 1 : numel(wtype)
       W = verticalWell(W, G, rock, wloc(1,w), wloc(2,w), 1, ...
           'Type', wtype{w}, 'Val', wtarget(w), ...
           'Radius', wrad(w), 'Name', wname{w}, ...
           'InnerProduct', 'ip_tpf', 'Comp_i', 1);
    end
    
    % --- Initiate and solve flow problem
    pv = poreVolume(G, rock);
    rS = initState(G, W, 0, 0.0);
    hT = computeTrans(G, rock);
    rS = incompTPFA(rS, G, hT, fluid, 'wells', W);
    for i=1:numel(rS.wellSol), rS.wellSol(i).sign = sgn(i); end
    
    % --- Compute RTD
    D   = computeTOFandTracer(rS, G, rock, 'wells', W, ...
        'maxTOF', inf, 'computeWellTOFs', true, 'firstArrival', true);
    WP  = computeWellPairs(rS, G, rock, W, D);
    rtd = computeRTD(rS, G, pv, D, WP, W, 'nsteps',2500);
    RTD = estimateRTD(pv, D, WP);
    
    % --- Plot RTD
    figure, hold on
    plot(rtd.t/year, rtd.values, '-', 'LineWidth', 2); set(gca,'ColorOrderIndex',1)
    plot(RTD.t/year, RTD.values, '--', 'LineWidth',2); set(gca,'ColorOrderIndex',1)
    tfa = min(D.ifa(W(2).cells))/year;      % first arrival
    tm  = rtd.volumes./rtd.allocations/year; % mean time
    vm  = max(RTD.values(:));
    plot([1; 1]*tm',repmat([0; vm],1,numel(rtd.volumes)), '-');
    legend('Simulated pulse: P1', 'Simulated pulse: P2', ...
        'From averaged TOF: P1', 'From averaged TOF: P2');    
    set(gca,'FontSize',14,'XLim',[0 130]); 
    xlabel('Time [years]')
    hold off
    
    axes('Position',[.6 .3 .3 .35]);
    plotCellData(G,log10(rock.perm(:,1)),'EdgeColor','none');
    outlineCoarseGrid(G,D.ppart,'LineWidth',1);
    view(2); axis equal off, colormap(parula)
    x = G.cells.centroids(vertcat(W.cells),:);
    hold on
    plot(x(:,1),x(:,2),'ok','MarkerSize',8,'MarkerFaceColor','w');
    text(x(:,1)+40,x(:,2), {W.name});
    hold off
    
    % --- F-Phi diagram
    [f, phi ] = computeFandPhi(rtd);
    [tf,tphi] = computeFandPhi(rtd, 'sum', true);

    figure, hold on
    plot(phi,f, '-', tphi, tf, '-','LineWidth',2);
    legend('Simulated pulse: P1', 'Simulated pulse: P2', ...
        'Simulated pulse: whole field','Location','SouthEast');
    hold off
    axis([0 1 0 1]); set(gca,'FontSize',14);
end