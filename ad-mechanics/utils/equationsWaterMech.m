function [eqs, names, types, state] = equationsWaterMech(state0, p, wellVars, state, model, dt, mechTerm, drivingForces, varargin)
%
%
% SYNOPSIS:
%   function [eqs, names, types, state] = equationsWaterMech(state0, p, wellVars, state, model, dt, mechTerm, drivingForces, varargin)
%
% DESCRIPTION: This function is very similar to equationsWater. The
% difference here is that it also takes as input mechanical terms, and the ADI
% initialization is not done here but by the model in the getEquations member
% function.
%
% PARAMETERS:
%   state0        - State at previous time step
%   p             - Pressure
%   wellVars      - Well variables
%   state         - State at given time step
%   model         - Model class instance that is used.
%   dt            - Time step
%   mechTerm      - Mechanical input which will enter the computation of the
%                   effective porosity
%   drivingForces - Structure that gathers the well parameters and boundary conditions.
%
% RETURNS:
%   eqs   - The residual values as ADI variables (that is with the Jacobian)
%           if the inputs were also ADI.
%   names - The name of each equations
%   types - The type of each equations
%   state - Some field related to well control of the state variables may be updated.
%
% EXAMPLE:
%
% SEE ALSO: equationsWater
%

    
    % Note that state is given only for output
    opt = struct('iteration', -1, ...
                 'resOnly', false); % just to avoid warning
    opt = merge_options(opt, varargin{:});

    % Shorter names for some commonly used parts of the model and forces.
    s = model.operators;
    f = model.fluid;
    rock = model.rock;
    G = model.G;
    W = drivingForces.W;

    [p0, wellSol0] = model.getProps(state0, 'pressure', 'wellSol');
    

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

    eqs{1} = (1 ./ dt) .*                                                            ...
             ((rock.poro .* (G.cells.volumes .* pvMult)  + rock.alpha .* mechTerm.new) .* bW -      ...
              (rock.poro .* (G.cells.volumes .* pvMult0) + rock.alpha .* mechTerm.old) .* f.bW(p0)) ...
             + s.Div(bWvW);

    names = {'water'};
    types = {'cell'};

    % Finally, add in and setup well equations
    wellSol = model.getProp(state, 'wellsol');
    [~, wellVarNames, wellMap] = ...
        model.FacilityModel.getAllPrimaryVariables(wellSol);
    [eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, ...
                                                      names, types, wellSol0, ...
                                                      wellSol, wellVars, ...
                                                      wellMap, p, {mobW}, ...
                                                      {rhoW}, {}, {}, dt, ...
                                                      opt);
    
end
