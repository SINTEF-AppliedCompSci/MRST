classdef OilWaterSurfactantModel1D < OilWaterSurfactantModel
   methods

      function model = OilWaterSurfactantModel1D(G, rock, fluid, varargin)
         model = model@OilWaterSurfactantModel(G, rock, fluid, varargin{:});
      end

      function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
         [state, report] = updateAfterConvergence@OilWaterSurfactantModel(model, state0, ...
                                                           state, dt, drivingForces);

         setFigure(1)
         clf
         subplot(2, 1, 1)
         plot(state.s(:, 1), '*-');
         axis([0, 100, 0, 1]);
         title('Saturation');
         subplot(2, 1, 2)
         plot(state.c, '*-');
         % axis([0, 100, 0, 50]);
         title('Concentration');
         drawnow

         setFigure(2)
         clf
         plot(state.c, '*-');
         axis([0, 100, 0, 51]);
         title('Concentration');
         drawnow

         setFigure(3)
         clf
         plot(state.ads, '*-');
         title('SURFADS');
         drawnow

         % setFigure(4)
         % clf
         % plot(state.SURFCNM, '*-');
         % title('SURFCNM');
         % drawnow

         % setFigure(5)
         % clf
         % plot(state.SURFST, '*-');
         % title('SURFST');
         % drawnow

         % waitforbuttonpress;
      end
   end
end
