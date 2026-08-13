classdef RockTransmissibilityTest < matlab.unittest.TestCase
    methods
        function test = RockTransmissibilityTest()
            mrstModule reset
            gravity('reset', 'off')
        end
    end

    methods (Test)
        function testMakeRockExpandsScalarAndNTG(test)
            G = computeGeometry(cartGrid([2, 2, 1], [2, 2, 1]));
            rock = makeRock(G, 100*milli*darcy, 0.25, 'ntg', 0.5);

            test.verifySize(rock.perm, [G.cells.num, 1]);
            test.verifySize(rock.poro, [G.cells.num, 1]);
            test.verifySize(rock.ntg, [G.cells.num, 1]);
            test.verifyEqual(unique(rock.poro), 0.25);
            test.verifyEqual(unique(rock.ntg), 0.5);
        end

        function testComputeTransInternalFaceSymmetry(test)
            G = computeGeometry(cartGrid([2, 1, 1], [2, 1, 1]));
            rock = makeRock(G, darcy, 0.2);
            T = computeTrans(G, rock);

            cf = G.cells.faces(:, 1);
            internalFace = cf(G.faces.neighbors(cf, 1) ~= 0 & G.faces.neighbors(cf, 2) ~= 0);
            internalFace = internalFace(1);
            halves = T(cf == internalFace);

            test.verifyEqual(numel(halves), 2);
            test.verifyGreaterThan(halves(1), 0);
            test.verifyEqual(halves(1), halves(2), 'AbsTol', 1e-12);
        end

        function testComputeTransAppliesNTGOnlyToLateralFaces(test)
            G = computeGeometry(cartGrid([2, 2, 2], [2, 2, 2]));
            rockBase = makeRock(G, darcy, 0.2);
            rockNTG  = makeRock(G, darcy, 0.2, 'ntg', 0.25);

            TBase = computeTrans(G, rockBase);
            TNTG  = computeTrans(G, rockNTG);
            dim = ceil(G.cells.faces(:, 2)./2);

            lateral = dim ~= 3;
            vertical = dim == 3;

            test.verifyEqual(TNTG(lateral), 0.25*TBase(lateral), 'RelTol', 1e-12);
            test.verifyEqual(TNTG(vertical), TBase(vertical), 'RelTol', 1e-12);
        end
    end
end
