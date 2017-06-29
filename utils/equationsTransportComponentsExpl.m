function [eqs, names, types] = equationsTransportComponentsExpl(state0, p, mastercomps, comps, state, ...
                                                      model, dt, drivingForces, ...
                                                      varargin)

    opt = struct('Verbose', mrstVerbose, ...
                 'reverseMode', false,...
                 'resOnly', false,...
                 'iteration', -1);  % Compatibility only

    opt = merge_options(opt, varargin{:});


    s = model.operators;
    G = model.G;
    f = model.fluid;

    chemModel = model.chemicalModel;
    nMC = chemModel.nMC;
    nC = chemModel.nC;
    mastercomps0 = cell(1, nMC);
    comps0 = cell(1, nC);
    [p0, mastercomps0{:}, comps0{:}] = model.getProps(state0, 'pressure', chemModel.MasterCompNames{:}, chemModel.CompNames{:});

    %Initialization of independent variables ----------------------------------

    %grav  = gravity;
    %gdz   = s.Grad(G.cells.centroids) * model.gravity';
    gdz   = s.Grad(G.cells.centroids) * model.getGravityVector()';

    %--------------------
    %check for p-dependent tran mult:
    trMult = 1;
    if isfield(f, 'tranMultR'), trMult = f.tranMultR(p); end

    %check for p-dependent porv mult:
    pvMult = 1; pvMult0 = 1;
    if isfield(f, 'pvMultR')
        pvMult =  f.pvMultR(p);
        pvMult0 = f.pvMultR(p0);
    end
    transMult=1;
    if isfield(f, 'transMult')
        transMult=f.transMult(p); 
    end

    trans=s.T.*transMult;
    % -------------------------------------------------------------------------
    % water props (calculated at oil pressure OK?)
    bW     = f.bW(p);
    rhoW   = bW.*f.rhoWS;
    
    % rhoW on face, avarge of neighboring cells (E100, not E300)
    rhoWf  = s.faceAvg(rhoW);
    mobW   = trMult./f.muW(p);
    dpW     = s.Grad(p) - rhoWf.*gdz;
    
    % water upstream-index
    upcw = (double(dpW)<=0);
    vW = - s.faceUpstr(upcw, mobW).*trans.*dpW;
    bWvW = s.faceUpstr(upcw, bW).*vW;

    if model.outputFluxes
        state = model.storeFluxes(state, vW, [], []);
    end

    if model.extraStateOutput
        state = model.storebfactors(state, bW, [], []);
        state = model.storeMobilities(state, mobW, [], []);
        state = model.storeUpstreamIndices(state, upcw, [], []);
    end

    % EQUATIONS ---------------------------------------------------------------
    % water:
    eqs{1} = (s.pv/dt).*( pvMult.*bW - pvMult0.*f.bW(p0) ) + s.Div(bWvW);

    chemModel  = model.chemicalModel;
    fluidMat   = model.fluidMat;
    nMC        = chemModel.nMC;
    surfMaster = chemModel.surfMaster;

    fluidConc = cell(1, nMC);
    
    for i = 1 : nMC
        if surfMaster(i)
            eqs{end + 1} = (s.pv/dt).*(mastercomps{i} - mastercomps0{i});
        else
            fluidConc{i} = 0;
            for j = 1 : nC
                fluidConc{i} = fluidConc{i} + fluidMat(i,j)*comps0{j};
            end
            eqs{end + 1} = (s.pv/dt).*(mastercomps{i}.*pvMult.*bW - mastercomps0{i}.* ...
                                       pvMult0.*f.bW(p0)) + s.Div(s.faceUpstr(upcw, fluidConc{i}).* bWvW);
        end
    end
                   
    names = {'water', chemModel.MasterCompNames{:}};

    types = cell(1, nMC + 1);
    [types{:}] = deal('cell');
    
    sW = ones(model.G.cells.num, 1); % dummy saturation
    [eqs, state, src] = model.addBoundaryConditionsAndSources(eqs, names, ...
                                                      types, state, {p}, {sW}, ...
                                                      {mobW}, {rhoW}, {}, ...
                                                      fluidConc, ...
                                                      drivingForces);
    
    
end

