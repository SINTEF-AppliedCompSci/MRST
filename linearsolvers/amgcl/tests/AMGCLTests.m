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

