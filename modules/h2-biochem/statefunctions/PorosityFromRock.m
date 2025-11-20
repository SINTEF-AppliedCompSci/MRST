classdef PorosityFromRock < StateFunction
    % Simple porosity provider for cases without bio-clogging.
    % Returns static porosity from model.rock.poro.

    methods
        function sf = PorosityFromRock(model, varargin)
            sf@StateFunction(model, varargin{:});
            sf.label = '\phi';
        end

        function phi = evaluateOnDomain(sf, model, state) %#ok<INUSD>
            phi = model.rock.poro;
        end
    end
end