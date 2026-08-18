classdef CoreADAndUtilitiesTest < matlab.unittest.TestCase
    methods
        function test = CoreADAndUtilitiesTest()
            mrstModule reset
        end
    end

    methods (Test, TestTags = {'tier1'})
        function initVariablesADISeedsIdentityJacobians(test)
            [x, y] = initVariablesADI(3, [1; 2]);

            test.verifyEqual(x.val, 3);
            test.verifyEqual(y.val, [1; 2]);
            test.verifyEqual(full(x.jac{1}), 1);
            test.verifyEqual(full(x.jac{2}), zeros(1, 2));
            test.verifyEqual(full(y.jac{1}), zeros(2, 1));
            test.verifyEqual(full(y.jac{2}), eye(2));
        end

        function adiArithmeticPropagatesDerivatives(test)
            [x, y] = initVariablesADI(2, 5);
            z = 3*x.*y + y.^2;

            test.verifyEqual(z.val, 55);
            test.verifyEqual(full(z.jac{1}), 15);
            test.verifyEqual(full(z.jac{2}), 16);
        end

        function simpleEquilibriumReturnsPartitionOfUnity(test)
            G = computeGeometry(cartGrid([1, 1, 4], [1, 1, 4]));
            s = simpleEquilibrium(G, 2);

            test.verifyEqual(size(s), [G.cells.num, 2]);
            test.verifyLessThan(max(abs(sum(s, 2) - 1)), 1e-12);
            test.verifyTrue(all(s(:) >= 0 & s(:) <= 1));
        end

        function simpleEquilibriumSupportsCustomVector(test)
            G = computeGeometry(cartGrid([4, 1, 1], [4, 1, 1]));
            s = simpleEquilibrium(G, 2, [1, 0, 0]);

            test.verifyEqual(size(s), [G.cells.num, 2]);
            test.verifyLessThan(max(abs(sum(s, 2) - 1)), 1e-12);
            test.verifyEqual(s(1, :), [0, 1], 'AbsTol', 1e-12);
            test.verifyEqual(s(end, :), [1, 0], 'AbsTol', 1e-12);
        end

        function simpleEquilibriumHandlesZeroWidthProjection(test)
            G = computeGeometry(cartGrid([1, 1, 1], [1, 1, 1]));
            G.faces.centroids(:, 1) = 0.5;
            G.cells.centroids(:, 1) = 0.5;

            s = simpleEquilibrium(G, 0.5, 1);

            test.verifyEqual(s, [0, 1]);
        end
    end

    methods (Test, TestTags = {'tier2'})
        function mrstModuleAndRequireEnforceDependencies(test)
            mrstModule clear
            test.verifyError(@() require('deckformat'), 'MRST:MissingModule');

            mrstModule add deckformat
            require('deckformat');

            mrstModule reset ad-core
            require('ad-core');
            test.verifyError(@() require('deckformat'), 'MRST:MissingModule');
        end
    end
end
