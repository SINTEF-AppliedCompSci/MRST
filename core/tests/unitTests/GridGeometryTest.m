classdef GridGeometryTest < matlab.unittest.TestCase
    methods
        function test = GridGeometryTest()
            mrstModule reset
        end
    end

    methods (Test)
        function testCartGridTopology2D(test)
            G = cartGrid([3, 2], [6, 4]);

            test.verifyEqual(G.griddim, 2);
            test.verifyEqual(G.cartDims(:).', [3, 2]);
            test.verifyEqual(G.cells.num, 6);
            test.verifyEqual(G.nodes.num, 12);
            test.verifyTrue(any(strcmp(G.type, 'tensorGrid')));
            test.verifyTrue(any(strcmp(G.type, 'cartGrid')));
        end

        function testComputeGeometryVolumesAndNormals(test)
            G = computeGeometry(cartGrid([2, 3, 4], [4, 6, 8]));

            test.verifyEqual(sum(G.cells.volumes), 4*6*8, 'AbsTol', 1e-10);
            test.verifyTrue(all(G.cells.volumes > 0));
            test.verifyEqual(sqrt(sum(G.faces.normals.^2, 2)), G.faces.areas, ...
                'AbsTol', 1e-10);
        end

        function testBoundaryFaceIndicesAliasesAndRanges(test)
            G = computeGeometry(cartGrid([4, 3, 2], [4, 3, 2]));

            leftFaces  = sort(boundaryFaceIndices(G, 'left'));
            westFaces  = sort(boundaryFaceIndices(G, 'west'));
            xminFaces  = sort(boundaryFaceIndices(G, 'xmin'));
            rangeFaces = sort(boundaryFaceIndices(G, 'left', 2:3, 1, 1:2));

            test.verifyEqual(leftFaces, westFaces);
            test.verifyEqual(leftFaces, xminFaces);
            test.verifyEqual(numel(leftFaces), 3*2);
            test.verifyEqual(numel(rangeFaces), 2);
            test.verifyTrue(all(ismember(rangeFaces, leftFaces)));
        end
    end
end
