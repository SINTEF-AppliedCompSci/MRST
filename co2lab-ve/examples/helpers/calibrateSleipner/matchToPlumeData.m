function obj = matchToPlumeData(model, states, plumes)

    numSteps = length(states); % in practice: either 1 or all the states

    obj = repmat({[]}, numSteps, 1);
    for step = 1:numSteps

        if numSteps == 1
            assert(~iscell(states)); % should just be a single state
            state = states; 
            plume = plumes;
        else
            state = states{step};
            plume = plumes{step};
        end
        
        if ~isempty(plume)
            sG = state.s(:,2);
            if iscell(sG)
                assert(length(sG) == 1)
                sG = sG{1};
            end
            o = matchToCo2Surface(sG, plumes{step}, model.G, model.fluid);
            alpha = 1/9; % tested alpha = 0, 0.1, 1 (a scaling value for dz)
            obj{step} = o + sum(alpha * (model.G.dz).^2 .* model.G.cells.volumes);
        else
            obj{step} = 0;
        end
    end
end
