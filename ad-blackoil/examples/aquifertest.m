%% load modules
mrstModule add mrst-gui deckformat ad-props ad-core ad-blackoil 

mrstVerbose true
gravity on

%% run two test cases

datadir = getDatasetPath('aquifertest', 'download', true);

for icase = 1 : 2
    
    fnroot = sprintf('2D_OILWATER_AQ_ECLIPSE_%d', icase);
    fn = fullfile(datadir, [fnroot, '.DATA']);
    [model, initState, schedule] = setupAquifertest(fn);
    [wellSols, states, schedulereport] = simulateScheduleAD(initState, model, ...
                                                      schedule);
    statesEcl = loadAquiferEclipseResult(datadir, fnroot);


    %% plots
    G = model.G;
    figure
    plotToolbar(G, states);
    view([2, 2]);
    colorbar
    title(sprintf('MRST computation - Case %d', icase));
    figure
    plotToolbar(G, statesEcl);
    view([2, 2]);
    colorbar
    title(sprintf('Eclipse computation - Case %d', icase));
    
    p = states{end}.pressure;
    pecl = statesEcl{end}.pressure;
    pecl = pecl*barsa;
    err = (p - pecl)./pecl;
    figure
    plotCellData(G, err);
    view([2, 2]);
    colorbar
    title(sprintf('Relative difference - Case %d', icase));
    
end
