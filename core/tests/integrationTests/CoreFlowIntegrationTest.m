classdef CoreFlowIntegrationTest < matlab.unittest.TestCase
    methods
        function test = CoreFlowIntegrationTest()
            mrstModule reset
            mrstModule add incomp
            gravity('reset', 'off')
        end
    end

    methods (Test)
        function testLinearPressureDropWithBoundaryConditions(test)
            G = computeGeometry(cartGrid([5, 1, 1], [5, 1, 1]));
            rock = makeRock(G, 100*milli*darcy, 1);
            fluid = initSingleFluid('mu', 1*centi*poise, ...
                                    'rho', 1000*kilogram/meter^3);
            T = computeTrans(G, rock);
            bc = [];
            bc = pside(bc, G, 'xmin', 2*barsa);
            bc = pside(bc, G, 'xmax', 1*barsa);

            state = incompTPFA(initResSol(G, 0), G, T, fluid, 'bc', bc);
            expected = interp1([0; 5], [2; 1]*barsa, G.cells.centroids(:, 1));

            test.verifyEqual(state.pressure, expected, 'RelTol', 1e-10);
            test.verifyTrue(all(diff(state.pressure) < 0));
        end

        function testGravityColumnPressureIsMonotone(test)
            gravity('reset', 'on')

            G = computeGeometry(cartGrid([1, 1, 12], [1, 1, 12]));
            rock = makeRock(G, 100*milli*darcy, 1);
            fluid = initSingleFluid('mu', 1*centi*poise, ...
                                    'rho', 1014*kilogram/meter^3);
            T = computeTrans(G, rock);
            bc = pside([], G, 'top', 100*barsa);

            state = incompTPFA(initResSol(G, 0), G, T, fluid, 'bc', bc);

            test.verifyGreaterThan(state.pressure(end), state.pressure(1));
            test.verifyTrue(all(diff(state.pressure) > 0));
            test.verifyEqual(state.facePressure(bc.face), 100*barsa, 'RelTol', 1e-12);
        end

        function testWellDrivenFlowProducesBalancedWellSolutions(test)
            G = computeGeometry(cartGrid([8, 1, 1], [8, 1, 1]));
            rock = makeRock(G, 100*milli*darcy, 0.25);
            fluid = initSingleFluid('mu', 1*centi*poise, ...
                                    'rho', 1000*kilogram/meter^3);
            T = computeTrans(G, rock);

            W = [];
            W = addWell(W, G, rock, 1, 'Type', 'rate', 'Val', 2*meter^3/day, ...
                'Sign', 1, 'WI', 1);
            W = addWell(W, G, rock, G.cells.num, 'Type', 'bhp', 'Val', barsa, ...
                'WI', 1);

            state = initState(G, W, barsa);
            state = incompTPFA(state, G, T, fluid, 'wells', W);

            q = cellfun(@sum, {state.wellSol.flux});
            test.verifyEqual(numel(state.wellSol), 2);
            test.verifyGreaterThan(q(1), 0);
            test.verifyLessThan(q(2), 0);
            test.verifyEqual(sum(q), 0, 'AbsTol', 1e-10);
        end
    end
end
