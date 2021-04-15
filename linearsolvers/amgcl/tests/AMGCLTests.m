classdef AMGCLTests < matlab.unittest.TestCase
    properties
        tolerance = 1e-12;
        checkAbsTol = 1e-8;
    end
    methods
        function test = AMGCLTests()
            mrstModule reset
            mrstModule add ad-unittest ad-core linearsolvers
        end
    end
    
    methods (Test)
        function scalarTest(test)
            [A, b, ref] = test.getSimpleMatrix(1);
            [x, err] = callAMGCL(A, b, 'tolerance', test.tolerance);
            test.assertEqual(x, ref, 'AbsTol', test.checkAbsTol)
            test.assertFalse(err > test.tolerance);
        end
        
        function multipleScalarRHSTest(test)
            [A, b, ref] = test.getSimpleMatrix();
            [x, err] = callAMGCL(A, b, 'tolerance', test.tolerance);
            test.assertEqual(x, ref, 'AbsTol', test.checkAbsTol)
            test.assertFalse(err > test.tolerance);
        end
        
        function blockTest(test)
            [A, b, ref] = test.getBlockMatrix(1);
            [x, err] = callAMGCL(A, b, 'tolerance', test.tolerance, 'block_size', 2);
            test.assertEqual(x, ref, 'AbsTol', test.checkAbsTol)
            test.assertFalse(err > test.tolerance);
        end
        
        function multipleBlockRHSTest(test)
            [A, b, ref] = test.getBlockMatrix(1);
            [x, err] = callAMGCL(A, b, 'tolerance', test.tolerance, 'block_size', 2);
            test.assertEqual(x, ref, 'AbsTol', test.checkAbsTol)
            test.assertFalse(err > test.tolerance);
        end
        
        function CPRScalarTest(test)
            [A, b, ref] = test.getBlockMatrix(1);
            block_size = 2;
            [x, err] = callAMGCL_cpr(A, b, block_size, 'cellMajorOrder', true, ....
                'cpr_blocksolver', false, ...
                'tolerance', test.tolerance, 'block_size', 2, 'use_drs', false);
            test.assertEqual(x, ref, 'AbsTol', test.checkAbsTol)
            test.assertFalse(err > test.tolerance);
        end
        
        function CPRBlockTest(test)
            [A, b, ref] = test.getBlockMatrix();
            block_size = 2;
            [x, err] = callAMGCL_cpr(A, b, block_size, 'cellMajorOrder', true, ....
                'cpr_blocksolver', true, ...
                'tolerance', test.tolerance, 'block_size', 2, 'use_drs', false);
            test.assertEqual(x, ref, 'AbsTol', test.checkAbsTol)
            test.assertFalse(err > test.tolerance);
        end
        
        function CPRDRSScalarTest(test)
            [A, b, ref] = test.getBlockMatrix(1);
            block_size = 2;
            [x, err] = callAMGCL_cpr(A, b, block_size, 'cellMajorOrder', true, ....
                'cpr_blocksolver', false, ...
                'tolerance', test.tolerance, 'block_size', 2, 'use_drs', true);
            test.assertEqual(x, ref, 'AbsTol', test.checkAbsTol)
            test.assertFalse(err > test.tolerance);
        end
        
        function CPRDRSBlockTest(test)
            [A, b, ref] = test.getBlockMatrix();
            block_size = 2;
            [x, err] = callAMGCL_cpr(A, b, block_size, 'cellMajorOrder', true, ....
                'cpr_blocksolver', true, ...
                'tolerance', test.tolerance, 'block_size', 2, 'use_drs', true);
            test.assertEqual(x, ref, 'AbsTol', test.checkAbsTol)
            test.assertFalse(err > test.tolerance);
        end
    end
    methods (Static)
        function [A, b, ref] = getSimpleMatrix(ix)
            A = sparse([1, 1, 1; 0, 1, 0; 0, 0, 1]);
            b = [1, 5, 12; 2, 3, 19; 3, 20, 19];
            if nargin > 0
                b = b(:, ix);
            end
            ref = A\b;
        end

        function [A, b, ref] = getBlockMatrix(ix)
            % Poisson-like block system
            A = sparse([ 3,  0, -1,  0, -1,  0;....
                         1, -1,  0,  0,  0,  0;
                        -1,  0,  2,  0, -1,  0;
                         0,  0,  1, -1,  0,  0;
                        -1,  0, -1,  0,  2,  0;
                         0,  0,  0,  0,  1, -1]);
            b = zeros(6, 1);
            b(1) = 10;
            b(end-1) = -10;
            n = 20;
            b = repmat(b, 1, n);
            b = bsxfun(@mtimes, b, (1:n));
            if nargin > 0
                b = b(:, ix);
            end
            ref = A\b;
        end
    end
end

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
