classdef SaturationDifferenceTransferFunction
    %SATURATION_DIFFERENCE_TRANSFER_FUNCTION
    
    properties
        beta;
        maximum_saturation;
        smooth_function_parameters;
    end
    
    methods
        function obj = SaturationDifferenceTransferFunction(beta,varargin)
            obj.maximum_saturation = 0.0;
            obj.smooth_function_parameters = [100,1];
            obj = merge_options(obj, varargin{:});
            obj.beta = beta;
        end
        
        function v = transfer(obj, model, state, domain_id)
            
            v = {};
            dom = model.G.FracturedDomains.domains{domain_id};
            
            %% Saturations
            s = model.getProps(state,'saturation');
                                    
            %% Matrix saturations and maximum matrix saturation 
            if obj.maximum_saturation > 0.0
                swmax = obj.maximum_saturation*ones(length(dom.connections),1);
                swf = s{1}(dom.connections(:,1));
                swm = s{1}(dom.connections(:,2));
            else
                swmax = s{1}(dom.connections(:,1));
                swm = s{1}(dom.connections(:,2));
            end
            
            %% Smooth function helps the numerical convergence
            % 'a' and 'b' are parameters controlling how fast the function
            % goes to 1. These are empirical parameters
            a = obj.smooth_function_parameters(1);
            b = obj.smooth_function_parameters(2);
            smooth_function = @(s) (1-exp(-sqrt(a/b).*s))/(1-exp(-sqrt(a/b)));
                        
            %% Fluxes
            vb = model.G.cells.volumes(dom.region);
            v{1} = vb.*obj.beta.*(swmax-swm);
            if obj.maximum_saturation > 0.0
                v{1} = vb.*obj.beta.*smooth_function(swf).*(swmax-swm);
            end
            v{2} = -v{1};          
        end
    end
end

