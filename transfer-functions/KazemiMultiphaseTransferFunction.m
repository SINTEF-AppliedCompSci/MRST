classdef KazemiMultiphaseTransferFunction
    %KAZEMI_MULTIPHASE_TRANSFER_FUNCTION
    
    properties
        shape_factor_model;
        gravity_model;
    end
    
    methods
        function obj = KazemiMultiphaseTransferFunction(shape_factor_model, varargin)
            obj.gravity_model = 'none';
            obj = merge_options(obj, varargin{:});
            obj.shape_factor_model = shape_factor_model;
        end
        
        function v = transfer(obj, model, state, domain_id)
            
            v = {};

            %% Get state variables
            [p, s, flag, mob, rho] = model.getProps(state,'PhasePressures',...
                                                          'Saturation',...
                                                          'PhaseUpwindFlag',...
                                                          'Mobility',...
                                                          'Density');
            
            %% Get domain object
            dom = model.G.FracturedDomains.domains{domain_id};
            
            %% Calculate shape factor
            sigma = obj.shape_factor_model.calculateShapeFactor(dom.rock);
            vb = model.G.cells.volumes(dom.region);
            
            %% computing term for gravity transfer
            if(~strcmp(obj.gravity_model,'none'))
                vg = {};
                lz = obj.shape_factor_model.block_dimension(:,3);
                swf = s{1}(dom.connections(:,1));
                swm = s{1}(dom.connections(:,2));
                
                hwf = swf.*lz; % TODO include residual saturations
                hwm = swm.*lz; % TODO include residual saturations 
                hof = lz - hwf;
                hom = lz - hwm;
                
                rhowf = rho{1}(dom.connections(:,1));
                rhoof = rho{2}(dom.connections(:,1));
                rhowm = rho{1}(dom.connections(:,2));
                rhoom = rho{2}(dom.connections(:,2));
                
                if(strcmp(obj.gravity_model,'eclipse_sonier'))
                    dwo = 0.5*9.81*(rhowm-rhoom);
                    vg{1} = - dwo.*(hwf-hwm);
                    vg{2} = + dwo.*(hwf-hwm);
                elseif(strcmp(obj.gravity_model,'gilman_kazemi'))
                    d = model.G.cells.centroids(dom.region,3);
                    vg{1} = - 9.81*rhowf.*(hwf-d) - 9.81*rhowm.*(hwm-d);
                    vg{2} = - 9.81*rhoof.*(hof-d) - 9.81*rhoom.*(hom-d);
                else
                    error(['Invalid Option in Gravity Model (',obj.gravity_model,')']);
                end
            end
            
            %% Compute fluxes (transfer)
            for i = 1:length(p)
                
                % Pressures in fracture and matrix
                pf{i} = p{i}(dom.connections(:,1));
                pm{i} = p{i}(dom.connections(:,2));
                
                % Fluxes
                fmob = model.operators.faceUpstr(flag{i}, mob{i});
                fmob = fmob(dom.global_connection_ids);
                v{i} = vb.*fmob.*sigma.*(pf{i}-pm{i});
                
                %% TODO: so far, it works only with oil/water systems
                if(model.water && model.oil && ~model.gas)
                    if(~strcmp(obj.gravity_model,'none'))
                        v{i} = v{i} + vb.*fmob.*sigma.*vg{i};
                    end
                end
            end
        end
    end
end

