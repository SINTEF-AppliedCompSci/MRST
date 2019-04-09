classdef ComponentPhaseMoleFractionsLV < GridProperty
    properties
    end
    
    methods
        function gp = ComponentPhaseMoleFractionsLV(model, varargin)
            gp@GridProperty(model, varargin{:});
        end

        function v = evaluateOnDomain(prop, model, state)
            [x, y] = model.getProps(state, 'x', 'y');
            if ~iscell(x)
                x = expandMatrixToCell(x);
            end
            if ~iscell(y)
                y = expandMatrixToCell(y);
            end
            if model.water
                ncomp = model.getNumberOfComponents();
                nc = model.G.cells.num;
                u = ones(nc, 1);
                x = [{[]}, x];
                y = [{[]}, y];
                w = cell(1, ncomp);
                w{1} = u;
                v = [w', x', y'];
            else
                v = [x', y'];
            end
        end
    end
end