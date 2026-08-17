classdef CoreRockFlowPrimitivesTest < matlab.unittest.TestCase
    methods
        function test = CoreRockFlowPrimitivesTest()
            mrstModule reset
        end
    end

    methods (Test, TestTags = {'tier1'})
        function makeRockExpandsScalarFields(test)
            G = computeGeometry(cartGrid([2, 1, 1], [2, 1, 1]));
            rock = makeRock(G, 100*milli*darcy, 0.25);
            pv   = poreVolume(G, rock);

            test.verifyEqual(size(rock.perm), [G.cells.num, 1]);
            test.verifyEqual(rock.poro, 0.25*ones(G.cells.num, 1));
            test.verifyEqual(pv, 0.25*ones(G.cells.num, 1), 'AbsTol', 1e-12);
        end

        function poreVolumeIncludesNetToGross(test)
            G = computeGeometry(cartGrid([2, 1, 1], [2, 1, 1]));
            rock = makeRock(G, [100; 200]*milli*darcy, [0.2; 0.4], 'ntg', [0.5; 1.0]);
            pv = poreVolume(G, rock);

            test.verifyEqual(pv, [0.1; 0.4], 'AbsTol', 1e-12);
        end

        function initResSolExpandsInputsAndAccountsForNNC(test)
            G = computeGeometry(cartGrid([2, 1, 1], [2, 1, 1]));
            G.nnc.cells = [1, 2];

            state = initResSol(G, 100, [0.2, 0.8]);

            test.verifyEqual(state.pressure, 100*ones(G.cells.num, 1));
            test.verifyEqual(state.s, repmat([0.2, 0.8], [G.cells.num, 1]));
            test.verifyEqual(size(state.flux), [G.faces.num + 1, 1]);
        end

        function initStateCreatesWellSolutions(test)
            G    = computeGeometry(cartGrid([2, 2, 2], [2, 2, 2]));
            rock = makeRock(G, 100*milli*darcy, 0.3);
            W    = addWell([], G, rock, 1, ...
                           'Type', 'bhp', ...
                           'Val', 100*barsa, ...
                           'Radius', 0.1, ...
                           'Dir', 'z', ...
                           'Name', 'W1', ...
                           'compi', [1, 0, 0]);

            state = initState(G, W, 100*barsa, [1, 0, 0]);

            test.verifyTrue(isfield(state, 'wellSol'));
            test.verifyEqual(numel(state.wellSol), 1);
            test.verifyEqual(state.wellSol(1).pressure, 100*barsa);
            test.verifyEqual(state.wellSol(1).flux, 0);
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

