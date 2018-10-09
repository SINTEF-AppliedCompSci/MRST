classdef SequentialPressureTransportModelDG < SequentialPressureTransportModel
    
    methods
        function model = SequentialPressureTransportModelDG(pressureModel, transportModel, varargin)
            model = model@SequentialPressureTransportModel(pressureModel, transportModel, varargin{:});
        end
        
        function [state, report] = stepFunction(model, state, state0, dt,...
                                                drivingForces, linsolve, nonlinsolve,...
                                                iteration, varargin)
                                            
            [state, report] = stepFunction@SequentialPressureTransportModel(model, state, state0, dt,...
                                                drivingForces, linsolve, nonlinsolve,...
                                                iteration, varargin{:});
            
            disc = model.transportModel.disc;
            
            if disc.limitAfterConvergence
                
                G = disc.G;
                transportModel = model.transportModel;
                [outside, jump] = deal(false);
                
                if disc.outTolerance < Inf
                    [sMin, sMax] = disc.getMinMaxSaturation(state);
                    outside = sMin < 0 - disc.outTolerance | ...
                              sMax > 1 + disc.outTolerance;
                end
                
                if disc.jumpTolerance < Inf
                    % Cells with interface jumps larger than threshold
                    [jumpVal, ~, cells] = disc.getInterfaceJumps(state.sdof(:,1), state);
                    j = accumarray(cells(:), repmat(jumpVal,2,1) > disc.jumpTolerance) > 0;
                    jump = false(G.cells.num,1);
                    jump(cells(:))          = j(cells(:));
                    jump(state.degree == 0) = false;
                end
                
                bad = outside | jump;
                
                if any(bad)
                    
                    [substate0, substate, submodel, subforces] = buildSubProblem(transportModel, state, state0, drivingForces, bad);
                    maps     = submodel.G.mappings.cellMap;
                    keepFull = maps.keep & ~submodel.G.parent.cells.ghost;
                    keepSub  = ~submodel.G.cells.ghost;
                    
                    % Make initial state equal to substate in ok cells and
                    % substate0 in trouble-cells
                    substate1 = substate;
                    ix0 = submodel.disc.getDofIx(substate0, Inf, keepSub, true);
                    ix  = submodel.disc.getDofIx(substate , Inf, keepSub, true);
                    substate1.sdof(ix(ix>0 & ix0>0),:) = substate0.sdof(ix0(ix>0 & ix0>0),:);
                    substate1.sdof(ix(ix>0 & ix0==0),:) = [];
                    substate1.degree(keepSub) = substate0.degree(keepSub);
                    substate1.s(keepSub,:) = substate0.s(keepSub,:);
                    substate.degree(keepSub) = 1;
                    
                    forceArg = getDrivingForces(submodel, subforces);
                    submodel.verbose = false;
                    nls = NonLinearSolver();
                    solveCell = find(submodel.G.mappings.cellMap.keep & ~submodel.G.parent.cells.ghost);
                    if numel(solveCell) == 1 || 1
                        fprintf('Solving cell %d ... ', solveCell);
                    end
                    [substate, ~] = ...
                        nls.solveTimestep(substate1, dt, submodel,...
                                    'initialGuess', substate, ...
                                    forceArg{:});

                    state.s(keepFull, :)   = substate.s(keepSub,:);
                    state.mob(keepFull, :) = substate.mob(keepSub,:);
                    state.degree(keepFull) = substate.degree(keepSub);
                    state.nDof(keepFull)   = substate.nDof(keepSub);

                    nDof = transportModel.disc.basis.nDof;
                    cells = (state.degree < transportModel.disc.degree) & keepFull;
                    ix = transportModel.disc.getDofIx(state, 2:nDof, cells);
                    state.sdof(ix,:) = [];    

                    state = transportModel.disc.updateDofPos(state);
                    state = transportModel.disc.mapDofs(state, state0);

                    ixFull = transportModel.disc.getDofIx(state, [], keepFull);
                    ixSub  = submodel.disc.getDofIx(substate, [], keepSub);
                    state.sdof(ixFull,:) = substate.sdof(ixSub,:);
                    
                end
            end
                
        end
        
    end
    
end