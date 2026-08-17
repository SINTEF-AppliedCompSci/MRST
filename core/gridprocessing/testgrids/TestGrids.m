classdef TestGrids < matlab.unittest.TestCase
    methods
        function test = TestGrids() %#ok<MANU>
            mrstModule reset
            mrstModule add incomp libgeometry triangle
        end
    end

    methods (Test)
        function cartesianGrdeclCanBeProcessed(test)
            xi = [0, 1, 3, 6];
            yi = [0, 2, 5];
            zi = [0, 1, 4];
            [X, Y] = ndgrid(linspace(0, 1, numel(xi)), ...
                            linspace(0, 1, numel(yi)));

            grdecl = cartesianGrdecl(xi, yi, zi, ...
                'depthz', sin(X).*cos(Y));

            G = test.processGrdecl(grdecl);

            test.verifyEqual(grdecl.cartDims, [3, 2, 2]);
            test.verifyEqual(G.cells.num, prod(grdecl.cartDims));
        end

        function simpleGrdeclSupportsNumericAndFunctionalFaults(test)
            dims = [4, 4, 3];

            numericGrdecl = simpleGrdecl(dims, 0, ...
                'flat', true, 'undisturbed', true, 'physDims', [2, 3, 4]);
            numericGrid = test.processGrdecl(numericGrdecl);

            test.verifyEqual(numericGrdecl.cartDims, dims);
            test.verifyEqual(max(numericGrid.nodes.coords), [2, 3, 4], ...
                'AbsTol', 1e-12);

            functionalGrdecl = simpleGrdecl(dims, @(x) 0.05 + 0.02*x);
            functionalGrid = test.processGrdecl(functionalGrdecl);

            test.verifyEqual(functionalGrdecl.cartDims, dims);
            test.verifyEqual(functionalGrid.cells.num, prod(dims));
        end

        function pinchedLayersGrdeclSupportsFaultVariants(test)
            rng(0);
            dims = [4, 3, 3];

            numericGrdecl = pinchedLayersGrdecl(dims, 0.1);
            numericGrid = test.processGrdecl(numericGrdecl);

            rng(0);
            functionalGrdecl = pinchedLayersGrdecl(dims, @(x) 0.02 + 0.03*x);
            functionalGrid = test.processGrdecl(functionalGrdecl);

            test.verifyEqual(numericGrdecl.cartDims, dims);
            test.verifyEqual(functionalGrdecl.cartDims, dims);
            test.verifyEqual(numericGrid.cells.num, prod(dims));
            test.verifyEqual(functionalGrid.cells.num, prod(dims));
        end

        function oneSlopingFaultRequiresEvenXDimension(test)
            didThrow = false;
            try
                oneSlopingFault([3, 2, 2], 0.1);
            catch ME
                didThrow = true;
                test.verifySubstring(ME.message, 'must be even');
            end
            test.verifyTrue(didThrow);

            grdecl = oneSlopingFault([4, 2, 2], 0.1);
            G = test.processGrdecl(grdecl);

            test.verifyEqual(grdecl.cartDims, [4, 2, 2]);
            test.verifyEqual(G.cells.num, prod(grdecl.cartDims));
        end

        function threeLayersBuildsLayeredPermeability(test)
            grdecl = threeLayers([2, 3], 3, [1, 2, 1]);
            G = test.processGrdecl(grdecl);

            test.verifyEqual(grdecl.cartDims, [5, 3, 4]);
            test.verifyEqual(unique(grdecl.PERMX), [1; 100]);
            test.verifyEqual(G.cells.num, prod(grdecl.cartDims));
        end

        function makeModel3BuildsProcessableFaultedModel(test)
            grdecl = makeModel3([6, 4, 3], [60, 40, 9]);
            G = test.processGrdecl(grdecl);

            test.verifyEqual(grdecl.cartDims, [6, 4, 3]);
            test.verifyGreaterThan(nnz(grdecl.ACTNUM == 0), 0);
            test.verifyEqual(G.cells.num, nnz(grdecl.ACTNUM));
        end

        function specialCaseCornerPointModelsCanBeProcessed(test)
            pinchedColumn = createPinchedColumn();
            pinchedColumnGrid = test.processGrdecl(pinchedColumn);
            test.verifyEqual(pinchedColumnGrid.cells.num, nnz(pinchedColumn.ACTNUM));

            raised = raisedColumn();
            raisedGrid = test.processGrdecl(raised);
            test.verifyEqual(raisedGrid.cells.num, prod(raised.cartDims));

            node = pinchedNode();
            nodeGrid = test.processGrdecl(node);
            test.verifyEqual(nodeGrid.cells.num, nnz(node.ACTNUM));
        end

        function pinchMiddleCellEnumeratesAllVariants(test)
            expectedCounts = [1, 4, 6, 4, 1];

            for k = 0:4
                grdecls = pinchMiddleCell(k);

                test.verifyNumElements(grdecls, expectedCounts(k + 1));
                for i = 1:numel(grdecls)
                    G = test.processGrdecl(grdecls{i});
                    test.verifyEqual(G.cells.num, prod(grdecls{i}.cartDims));
                end
            end
        end

        function twisterSupportsPointsAndGridStructures(test)
            pt = [0, 0; 1, 0; 0, 1; 1, 1];
            twistedPt = twister(pt, 0.01);
            test.verifySize(twistedPt, size(pt));
            test.verifyNotEqual(twistedPt, pt);

            G = cartGrid([3, 3, 2], [1, 1, 1]);
            coords = G.nodes.coords;
            G = twister(G, 0.01);
            test.verifySize(G.nodes.coords, size(coords));
            test.verifyNotEqual(G.nodes.coords, coords);

            G = computeGeometry(G);
            test.verifyTrue(checkGrid(G));
        end

        function extrudedTriangleGridSupportsPrismaticAndDualGrids(test)
            G1 = extrudedTriangleGrid(500);
            G1 = computeGeometry(G1);
            test.verifyTrue(checkGrid(G1));
            test.verifyEqual(G1.griddim, 3);
            test.verifyGreaterThan(G1.cells.num, 0);

            G2 = extrudedTriangleGrid(500, true);
            G2 = computeGeometry(G2);
            test.verifyTrue(checkGrid(G2));
            test.verifyEqual(G2.griddim, 3);
            test.verifyGreaterThan(G2.cells.num, 0);
        end

        function cutCellGridWorkflowProducesValidCutAndSliceGrids(test)
            G = addBoundingBoxFields(computeGeometry(cartGrid([4, 4, 2], [4, 4, 2])));

            [G, ix, g2] = sliceGrid(G, [2, 2, 1], ...
                'normal', [1, -1, 0.5], 'radius', inf, 'topoSplit', true);

            test.verifyNotEmpty(ix);
            test.verifyGreaterThan(nnz(ix.new.faces == 3), 0);
            test.verifyTrue(checkGrid(G));

            g2 = repairNormals(computeGeometry(g2));
            test.verifyTrue(checkGrid(g2));
            test.verifyEqual(g2.griddim, 2);

            points = [0, 2, 0.2; 1, 2.5, 0.6; 2, 2, 1.0; 4, 1.5, 1.6];
            [Gcurve, ixCurve, gcurve] = sliceGrid( ...
                addBoundingBoxFields(computeGeometry(cartGrid([6, 6, 3], [6, 6, 3]))), ...
                points, 'cutDir', [0, -0.5, 1], 'radius', inf);

            test.verifyGreaterThan(nnz(ixCurve.new.faces == 3), 0);
            test.verifyTrue(checkGrid(Gcurve));

            gcurve = repairNormals(computeGeometry(gcurve));
            test.verifyTrue(checkGrid(gcurve));
            test.verifyEqual(gcurve.griddim, 2);
        end

        function markCutGridsExampleOneProducesExpectedRegions(test)
            nx = 5;
            physdims = [1, 1, 1];
            G0 = computeGeometry(cartGrid([nx, nx, nx], physdims));
            n = [1, 1, 0];

            [Gsingle, gixSingle] = sliceGrid(G0, physdims/2, 'normal', n);
            singleDomains = markCutGrids(Gsingle, gixSingle.new.faces);
            test.verifyEqual(numel(unique(singleDomains)), 2);

            cuts = [physdims/2; physdims/2 + 0.1];
            [Gdouble, gixDouble] = sliceGrid(G0, cuts, 'normal', n);
            doubleDomains = markCutGrids(Gdouble, gixDouble.new.faces);
            test.verifyEqual(numel(unique(doubleDomains)), 3);
        end

        function markCutGridsExampleTwoTracksMappedFaces(test)
            rng(0);
            G = addBoundingBoxFields(computeGeometry(cartGrid([5, 5, 5], [20, 20, 20])));
            nCuts = 3;
            parentFaces = cell(nCuts, 1);

            for k = 1:nCuts
                [G, gix] = sliceGrid(G, 10 + 6*rand(1, 3) - 3, ...
                    'normal', 2*rand(1, 3) - 1, ...
                    'radius', inf, 'topoSplit', true);
                test.verifyTrue(checkGrid(G));
                parentFaces{k} = gix.parent.faces;
            end

            faceStatus = gix.new.faces;
            for m = 2:nCuts
                fk = find(parentFaces{m - 1} == 0);
                for k = m:nCuts
                    fk = test.mapFace(parentFaces{k}, fk);
                end
                faceStatus(fk) = 3;
            end

            domains = markCutGrids(G, faceStatus);
            domainIds = unique(domains);

            test.verifyGreaterThan(numel(domainIds), 1);
            test.verifyEqual(domainIds(:), (1:numel(domainIds)).');
        end
    end

    methods (Access = private)
        function G = processGrdecl(test, grdecl)
            G = computeGeometry(processGRDECL(grdecl));
            test.verifyTrue(checkGrid(G));
        end

        function fmap = mapFace(~, parentFaces, f)
            fmap = [];
            for i = 1:numel(f)
                fmap = [fmap; find(parentFaces == f(i))]; %#ok<AGROW>
            end
        end
    end
end

%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
