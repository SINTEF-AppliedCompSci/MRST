classdef BoundaryAndStateTest < matlab.unittest.TestCase
    methods
        function test = BoundaryAndStateTest()
            mrstModule reset
            gravity('reset', 'off')
        end
    end

    methods (Test)
        function testPsideSelectsExpectedFaces(test)
            G = computeGeometry(cartGrid([3, 2, 1], [3, 2, 1]));
            faces = boundaryFaceIndices(G, 'left');
            bc = pside([], G, 'left', 100*barsa, 'sat', [1, 0]);

            test.verifyEqual(sort(bc.face), sort(faces));
            test.verifyEqual(numel(bc.type), numel(faces));
            test.verifyTrue(all(strcmp(bc.type, 'pressure')));
            test.verifyEqual(bc.value, repmat(100*barsa, numel(faces), 1));
            test.verifyEqual(size(bc.sat), [numel(faces), 2]);
        end

        function testFluxsideSplitsTotalFluxByArea(test)
            G = computeGeometry(cartGrid([1, 2, 1], [1, 2, 1]));
            totalFlux = 4*meter^3/day;
            bc = fluxside([], G, 'left', totalFlux, 'sat', 1);

            test.verifyEqual(sum(bc.value), totalFlux, 'AbsTol', 1e-12);
            test.verifyEqual(numel(bc.face), 2);
            test.verifyEqual(bc.value(1), bc.value(2), 'AbsTol', 1e-12);
        end

        function testAddWellLogicalMaskAndInitState(test)
            G = computeGeometry(cartGrid([3, 1, 1], [3, 1, 1]));
            rock = makeRock(G, 100*milli*darcy, 0.3);
            mask = [true; false; true];

            W = addWell([], G, rock, mask, ...
                'Type', 'rate', 'Val', meter^3/day, ...
                'WI', ones(nnz(mask), 1));
            state = initState(G, W, 50*barsa, [1, 0, 0]);

            test.verifyEqual(W.cells, find(mask));
            test.verifyEqual(numel(state.wellSol), 1);
            test.verifyEqual(state.wellSol(1).pressure, 50*barsa);
            test.verifyEqual(state.wellSol(1).flux, zeros(nnz(mask), 1));
            test.verifyEqual(size(state.s), [G.cells.num, 3]);
        end

        function testInitStateRejectsPhaseMismatch(test)
            G = computeGeometry(cartGrid([2, 1, 1], [2, 1, 1]));
            rock = makeRock(G, 100*milli*darcy, 0.3);
            W = addWell([], G, rock, 1, ...
                'Type', 'rate', 'Val', meter^3/day, ...
                'WI', 1, 'compi', [1, 0, 0]);

            test.verifyError(@() initState(G, W, barsa, [1, 0]), ...
                'MATLAB:assertion:failed');
        end
    end
end
