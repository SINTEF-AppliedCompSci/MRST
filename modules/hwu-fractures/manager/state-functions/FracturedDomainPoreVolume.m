classdef FracturedDomainPoreVolume < StateFunction
    % Effective pore-volume after pressure-multiplier
    properties
    end
    
    methods
        function gp = FracturedDomainPoreVolume(model, varargin)
            gp@StateFunction(model, varargin{:});
            if isfield(model.fluid, 'pvMultR')
                gp = gp.dependsOn({'pressure'}, 'state');
            end
        end
        function pv = evaluateOnDomain(prop, model, state)
            % Get effective pore-volume, accounting for possible
            % rock-compressibility
            f = model.fluid;
            pv = model.operators.pv;
            if isfield(f, 'pvMultR')
                p = model.getProp(state, 'pressure');
                pvMult = prop.evaluateFunctionOnDomainWithArguments(f.pvMultR, p);
                pv = pv.*pvMult;
            end
            for i = 1:length(model.G.FracturedDomains.domains)
                
                dom = model.G.FracturedDomains.domains{i};
                if strcmp(dom.type,'multi_continuum')
                    v = model.G.cells.volumes(dom.region);
                    if isfield(dom.fluid, 'pvMultR')
                        p = model.getProp(state, 'pressure');
                        pvMult = prop.evaluateFunctionOnDomainWithArguments(dom.fluid.pvMultR, p);
                        v = v.*pvMult(dom.connection_cell_list);                    
                    end                    
                else
                    ap = dom.rock.aperture;
                    area = model.G.faces.areas(dom.region);
                    v = ap.*area;
                end
                pvdom = dom.rock.poro.*v;
                pv = [pv; pvdom];
            end
        end
    end
end