classdef CoreGridGeometryTest < matlab.unittest.TestCase
    methods
        function test = CoreGridGeometryTest()
            mrstModule reset
        end
    end

    methods (Test, TestTags = {'tier1'})
        function cartGridConstructsExpectedTopology(test)
            G = cartGrid([2, 3, 1], [4, 6, 2]);

            test.verifyEqual(G.cartDims, [2, 3, 1]);
            test.verifyEqual(G.griddim, 3);
            test.verifyEqual(G.cells.num, 6);
            test.verifyEqual(G.faces.num, 29);
            test.verifyEqual(G.nodes.num, 24);
            test.verifyEqual(diff(G.cells.facePos), repmat(6, [G.cells.num, 1]));
            test.verifyEqual(size(G.faces.neighbors), [G.faces.num, 2]);
            test.verifyTrue(any(any(G.faces.neighbors == 0, 2)));
        end

        function tensorGridAndRemoveCellsMaintainConnectivity(test)
            G = tensorGrid([0, 1, 3], [0, 2, 5]);
            G = computeGeometry(G);

            [H, cellmap, facemap, nodemap] = removeCells(G, 2);

            test.verifyEqual(H.cells.num, 3);
            test.verifyEqual(cellmap, [1; 3; 4]);
            test.verifyEqual(numel(facemap), H.faces.num);
            test.verifyEqual(numel(nodemap), H.nodes.num);
            test.verifyTrue(all(H.cells.faces(:, 1) >= 1 & H.cells.faces(:, 1) <= H.faces.num));
            test.verifyTrue(all(H.faces.nodes >= 1 & H.faces.nodes <= H.nodes.num));
        end

        function triangleGridGeometryIsConsistent(test)
            p = [0, 0; ...
                 1, 0; ...
                 1, 1; ...
                 0, 1];
            t = [1, 2, 3; ...
                 1, 3, 4];

            G = triangleGrid(p, t);
            G = computeGeometry(G);

            test.verifyEqual(G.cells.num, 2);
            test.verifyEqual(G.faces.num, 5);
            test.verifyEqual(sum(G.cells.volumes), 1, 'AbsTol', 1e-12);
            test.verifyEqual(size(G.cells.centroids), [2, 2]);
            CoreGridGeometryTest.verifyFaceNormalLengthsMatchAreas(G, 1e-12);
        end

        function tetrahedralGridGeometryIsConsistent(test)
            p = [0, 0, 0; ...
                 1, 0, 0; ...
                 0, 1, 0; ...
                 0, 0, 1];
            t = [1, 2, 3, 4];

            G = tetrahedralGrid(p, t);
            G = computeGeometry(G);

            test.verifyEqual(G.cells.num, 1);
            test.verifyEqual(G.faces.num, 4);
            test.verifyEqual(G.cells.centroids, [0.25, 0.25, 0.25], 'AbsTol', 1e-12);
            test.verifyEqual(G.cells.volumes, 1/6, 'AbsTol', 1e-12);
            CoreGridGeometryTest.verifyFaceNormalLengthsMatchAreas(G, 1e-12);
        end

        function computeGeometryUpdatesAfterNodeChanges(test)
            G = cartGrid([2, 2], [2, 2]);
            G = computeGeometry(G);
            centroids0 = G.cells.centroids;

            G.nodes.coords(:, 1) = G.nodes.coords(:, 1) + 0.25*G.nodes.coords(:, 2);
            G = computeGeometry(G);

            test.verifyGreaterThan(norm(G.cells.centroids - centroids0, 'fro'), 0);
            CoreGridGeometryTest.verifyFaceNormalLengthsMatchAreas(G, 1e-12);
        end
    end

    methods (Static, Access = private)
        function verifyFaceNormalLengthsMatchAreas(G, tol)
            n = sqrt(sum(G.faces.normals.^2, 2));
            assert(all(abs(n - G.faces.areas) <= tol));
        end
    end
end
