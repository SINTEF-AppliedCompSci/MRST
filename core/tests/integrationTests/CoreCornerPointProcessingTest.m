classdef CoreCornerPointProcessingTest < matlab.unittest.TestCase
    methods
        function test = CoreCornerPointProcessingTest()
            mrstModule reset
        end
    end

    methods (Test, TestTags = {'tier2'})
        function processFaultedSimpleGrdeclProducesTaggedFaces(test)
            grdecl = simpleGrdecl([4, 3, 2], 0.15, 'flat', true);
            G      = processGRDECL(grdecl);
            G      = computeGeometry(G);

            test.verifyEqual(numel(G), 1);
            test.verifyEqual(G.cells.num, prod(grdecl.cartDims));
            test.verifyTrue(any(G.faces.tag > 0));
            test.verifyTrue(all(G.cells.volumes > 0));
            CoreCornerPointProcessingTest.verifyFaceNormalLengthsMatchAreas(G, 1e-10);
        end

        function processGrdeclRespectsActnum(test)
            grdecl = simpleGrdecl([3, 2, 2], 0, 'flat', true, 'undisturbed', true);
            actnum = true(prod(grdecl.cartDims), 1);
            actnum([2, 5]) = false;
            grdecl.ACTNUM = int32(actnum);

            G = processGRDECL(grdecl, 'SplitDisconnected', false);
            G = computeGeometry(G);

            test.verifyEqual(G.cells.num, nnz(actnum));
            test.verifyEqual(sort(G.cells.indexMap), find(actnum));
            test.verifyTrue(all(G.cells.volumes > 0));
        end

        function splitDisconnectedCornerPointGridReturnsComponents(test)
            grdecl = simpleGrdecl([3, 1, 1], 0, 'flat', true, 'undisturbed', true);
            grdecl.ACTNUM = int32([1; 0; 1]);

            G = processGRDECL(grdecl, 'CheckGrid', false);

            test.verifyEqual(numel(G), 2);
            test.verifyEqual(sort([G(1).cells.num, G(2).cells.num]), [1, 1]);
        end

        function processedThreeLayerGridSupportsDownstreamTransmissibility(test)
            grdecl = threeLayers([1, 2], 2, [1, 2, 1]);
            G      = processGRDECL(grdecl, 'SplitDisconnected', false);
            G      = computeGeometry(G);
            rock   = makeRock(G, grdecl.PERMX(G.cells.indexMap)*milli*darcy, 0.2);
            hT     = computeTrans(G, rock);

            test.verifyEqual(size(hT), [size(G.cells.faces, 1), 1]);
            test.verifyTrue(all(isfinite(hT)));
            test.verifyTrue(all(hT > 0));
        end
    end

    methods (Static, Access = private)
        function verifyFaceNormalLengthsMatchAreas(G, tol)
            n = sqrt(sum(G.faces.normals.^2, 2));
            assert(all(abs(n - G.faces.areas) <= tol));
        end
    end
end
