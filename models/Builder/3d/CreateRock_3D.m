function model = CreateRock_3D(model)
    experiment = model.experiment;
    simulation = model.simulation;
    grid       = model.grid;
    G          = grid.G;
    poro = experiment.rock.poro.value * ones(G.cells.num, 1);
    perm = experiment.rock.perm.value * ones(G.cells.num, 1); 

%     poro = bsxfun(@plus,poro,randn(numel(poro),1).*poro*0.5);
%     poro(poro < 0.001) = 0.001; poro(poro > 1) = 1;
    perm = bsxfun(@plus,perm,randn(numel(perm),1).*perm*1);
    
%     poro = use_porosity_profile(simulation, G, experiment);
%     figure; plot(G.cells.centroids(:,1), poro);
       
    rock = makeRock(G, perm, poro); 
    rock.pv = poreVolume(G, rock); 
    if (isfield(simulation,'bCells'))
        rock.regions = struct('saturation', grid.satNum);
    end
    model.rock = rock;
end