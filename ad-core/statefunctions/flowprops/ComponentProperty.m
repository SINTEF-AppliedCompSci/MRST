classdef ComponentProperty
    properties
        
    end
    
    methods
        function gp = ComponentProperty(model)
            if nargin > 0
                ncomp = model.getNumberOfComponents();
                deps = cell(ncomp, 1);
                exts = cell(ncomp, 1);
                for c = 1:ncomp
                    deps{c} = model.Components{c}.dependencies;
                    exts{c} = model.Components{c}.externals;
                end
                % Internal - flow props dependencies
                deps = unique(vertcat(deps{:}));
                gp = gp.dependsOn(deps); %#ok virtual class
                % Manage external dependencies
                exts = vertcat(exts{~cellfun(@isempty, exts)});
                if ~isempty(exts)
                    names = {exts.name};
                    [~, pos] = unique(names);
                    exts = exts(pos);
                    gp = gp.dependsOn(exts); %#ok virtual class
                end
            end
        end
    end
end