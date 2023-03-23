classdef Strain < StateFunction
    % Compute strain, depending on arguments make this cell-wise or over
    % the domain
    properties
    end
    
    methods
        function gp = Strain(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'xd'}, 'state');
        end
        function v = evaluateOnDomain(prop, model, state)
            mechModel = model.mechModel;
            u = mechModel.operators.V_dir; 
            u(~mechModel.operators.isdirdofs) = model.getProp(state, 'xd');
            
            if(model.G.griddim == 2)
                lin_dim = 3; 
            elseif(model.G.griddim == 3)
                lin_dim = 6;
            else
                error('Wrong dimsension');
            end
            e = reshape(model.operators.global_strain * u, lin_dim, [])';
            
            try 
                switch mechModel.mech.deformation_reference
                    case 'cells'
                        v = e;
                    case 'domain'
                        v = e.*(model.G.cells.volumes...
                                ./prod(max(Dmodel.G.nodes.coords)));
                end
            catch
                % if nothing is specified just return cell-wise strains
                v = e;
            end
        end
    end
end