classdef OilWaterSurfactantModel1D < FullyImplicitOilWaterSurfactantModel

    methods
        function model = OilWaterSurfactantModel1D(G, rock, fluid, varargin)
            model = model@FullyImplicitOilWaterSurfactantModel(G, rock, fluid, varargin{:});
        end

        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@FullyImplicitOilWaterSurfactantModel(model, state0, ...
                                                              state, dt, drivingForces);

            setFigure(1)
            clf
            subplot(2, 1, 1)
            plot(state.s(:, 1), '-');
            title('Saturation');
            subplot(2, 1, 2)
            plot(state.c, '-');
            % axis([0, 100, 0, 50]);
            title('Concentration');
            drawnow

        end
    end
end
