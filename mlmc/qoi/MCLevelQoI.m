classdef MCLevelQoI < BaseQoI
    
    properties
        levelQoIs
    end
    
    methods
        function qoi = MCLevelQoI(levelQoIs)
            cls = cellfun(@class, levelQoIs, 'UniformOutput', false);
            assert(all(strcmpi(cls{1}, cls)));
            assert(numel(levelQoIs) <= 2);
            qoi = qoi@BaseQoI();
            qoi.levelQoIs = levelQoIs;
            qoi.names = qoi.levelQoIs{1}.names;
        end
        
        function du = getQoI(qoi, seed)
            du = qoi.levelQoIs{end}.ResultHandler{seed};
            if numel(qoi.levelQoIs) == 2
                u0 = qoi.levelQoIs{1}.ResultHandler{seed};
                for fn = qoi.names
                    du.(fn{1}) = du.(fn{1}) - u0.(fn{1});
                end
                du.cost = du.cost + u0.cost;
            end
            qoi.ResultHandler{seed} = du;
        end
        
    end
end