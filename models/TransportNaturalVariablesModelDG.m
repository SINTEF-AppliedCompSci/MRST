classdef TransportNaturalVariablesModelDG < TransportNaturalVariablesModel
    properties
        disc
        tryMaxDegree
    end
    
    methods
        function model = TransportNaturalVariablesModelDG(G, rock, fluid, compFluid, varargin)
            
            model = model@TransportNaturalVariablesModel(G, rock, fluid, compFluid, varargin{:});
            model.tryMaxDegree = true;
        end

        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = transportEquationsNaturalVarsDG(state0, state, model, dt, ...
                            drivingForces, varargin{:});
        end
        
        % ----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name)
            % Map variables to state field.
            %
            % SEE ALSO:
            %   :meth:`ad_core.models.PhysicalModel.getVariableField`
            switch(lower(name))
                case {'swdof'}
                    index = 1;
                    fn = 'sdof';
                case {'sodof'}
                    index = 2;
                    fn = 'sdof';
                case {'sgdof'}
                    index = 3;
                    fn = 'sdof';
                case {'sdof'}
                    index = ':';
                    fn = 'sdof';
                case {'xdof'}
                    index = ':';
                    fn = 'xdof';
                case {'ydof'}
                    index = ':';
                    fn = 'ydof';
                case {'componentsdof'}
                    index = ':';
                    fn = 'componentsdof';
                otherwise
                    % This will throw an error for us
                    [fn, index] = getVariableField@TransportNaturalVariablesModel(model, name);
            end
        end
        
        %{
        % ----------------------------------------------------------------%
        function [xM,  yM,  rhoO,  rhoG,  muO,  muG, f_L, f_V,...
                  xM0, yM0, rhoO0, rhoG0] = ...
                  getTimestepPropertiesEoS(model, state, state0, p, temp, xdof, ydof, zdof, sOdof, sGdof, cellJacMap)
            
            eos = model.EOSModel;
            ncomp = numel(xdof);
%             % Get cell averages of dg variables
%             [x, y] = deal(cell(ncomp,1));
%             z = zeros(model.G.cells.num, ncomp);
%             for cNo = 1:ncomp
%                 x{cNo}   = model.disc.getCellMean(xdof{cNo}, state);
%                 y{cNo}   = model.disc.getCellMean(ydof{cNo}, state);
%                 z(:,cNo) = model.disc.getCellMean(zdof(:,cNo), state);
%             end
%             sO = model.disc.getCellMean(sOdof, state);
%             sG = model.disc.getCellMean(sGdof, state);
            
            [Z_L, Z_V, f_L, f_V] = eos.getProperties(p, temp, x, y, z, sO, sG, state, cellJacMap);

            [xM, rhoO, muO] = model.getFlowPropsNatural(p, x, Z_L, temp, true);
            [yM, rhoG, muG] = model.getFlowPropsNatural(p, y, Z_V, temp, false);
            
            [p0, temp0, x0, y0] = model.getProps(state0, 'pressure', 'T', 'x', 'y');
            x0 = expandMatrixToCell(x0);
            y0 = expandMatrixToCell(y0);

            [xM0, rhoO0] = model.getFlowPropsNatural(p0, x0, state0.Z_L, temp0, true);
            [yM0, rhoG0] = model.getFlowPropsNatural(p0, y0, state0.Z_V, temp0, false);

        end
        %}
        
        % ----------------------------------------------------------------%
        
        
    end
end