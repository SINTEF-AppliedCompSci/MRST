classdef CoreRockFlowPrimitivesIntegrationTest < matlab.unittest.TestCase
    methods
        function test = CoreRockFlowPrimitivesIntegrationTest()
            mrstModule reset
        end
    end

    methods (Test, TestTags = {'tier2'})
        function structuredFlowSetupProducesPositiveTransmissibilities(test)
            G    = computeGeometry(cartGrid([2, 1, 1], [2, 1, 1]));
            rock = makeRock(G, 1, 0.2);
            hT   = computeTrans(G, rock);

            cf = G.cells.faces(:, 1);
            nf = G.faces.num;
            T  = 1 ./ accumarray(cf, 1 ./ hT, [nf, 1]);
            intInx = all(G.faces.neighbors ~= 0, 2);

            test.verifyEqual(size(hT), [size(G.cells.faces, 1), 1]);
            test.verifyTrue(all(isfinite(hT)));
            test.verifyTrue(all(hT > 0));
            test.verifyEqual(nnz(intInx), 1);
            test.verifyEqual(T(intInx), 1, 'AbsTol', 1e-12);
        end
    end
end
