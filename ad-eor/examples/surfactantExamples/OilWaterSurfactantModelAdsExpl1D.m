classdef OilWaterSurfactantModelAdsExpl1D < OilWaterSurfactantModelAdsExpl
    methods

        function model = OilWaterSurfactantModelAdsExpl1D(G, rock, fluid, varargin)
            model = model@OilWaterSurfactantModelAdsExpl(G, rock, fluid, varargin{:});
        end

        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@OilWaterSurfactantModelAdsExpl(model, state0, ...
                                                              state, dt, drivingForces);
            try
                set(0, 'CurrentFigure', 1);
            catch
                figure(1);
            end
            subplot(2, 1, 1)
            plot(state.s(:, 1));
            axis([0, 100, 0, 1]);
            title('Saturation');
            subplot(2, 1, 2)
            plot(state.c);
            axis([0, 100, 0, 50]);
            title('Concentration');
            drawnow

        end
    end
end
