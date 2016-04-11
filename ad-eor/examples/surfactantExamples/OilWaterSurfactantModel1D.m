classdef OilWaterSurfactantModel1D < ImplicitExplicitOilWaterSurfactantModel
    methods

        function model = OilWaterSurfactantModel1D(G, rock, fluid, varargin)
            model = model@ImplicitExplicitOilWaterSurfactantModel(G, rock, fluid, varargin{:});
        end

        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@ImplicitExplicitOilWaterSurfactantModel(model, state0, ...
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

            setFigure(2)
            clf
            plot(state.c, '-');
            title('Concentration');
            drawnow

            setFigure(3)
            clf
            plot(state.ads, '-');
            title('SURFADS');
            drawnow

        end
    end
end
