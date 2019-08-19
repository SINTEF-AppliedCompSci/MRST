classdef AMGCLTests < matlab.unittest.TestCase
    methods
        function test = AMGCLTests()
            mrstModule reset
            mrstModule add ad-unittest ad-core linearsolvers
        end
    end
    
    methods (Test)
        function scalarTest(test)
            [A, b, ref] = test.getSimpleMatrix(1);
            tol = 1e-8;
            [x, err] = callAMGCL(A, b, 'tolerance', tol);
            test.assertEqual(x, ref, 'AbsTol', 10*tol)
            test.assertFalse(err > tol);
        end
        
        function multipleScalarRHSTest(test)
            [A, b, ref] = test.getSimpleMatrix();
            tol = 1e-8;
            [x, err] = callAMGCL(A, b, 'tolerance', tol);
            test.assertEqual(x, ref, 'AbsTol', 10*tol)
            test.assertFalse(err > tol);
        end
        
        function blockTest(test)
            [A, b, ref] = test.getBlockMatrix(1);
            tol = 1e-8;
            [x, err] = callAMGCL(A, b, 'tolerance', tol, 'block_size', 2);
            test.assertEqual(x, ref, 'AbsTol', 10*tol)
            test.assertFalse(err > tol);
        end
        
        function multipleBlockRHSTest(test)
            [A, b, ref] = test.getBlockMatrix(1);
            tol = 1e-8;
            [x, err] = callAMGCL(A, b, 'tolerance', tol, 'block_size', 2);
            test.assertEqual(x, ref, 'AbsTol', 10*tol)
            test.assertFalse(err > tol);
        end
        
        function CPRScalarTest(test)
            [A, b, ref] = test.getBlockMatrix(1);
            tol = 1e-8;
            block_size = 2;
            [x, err] = callAMGCL_cpr(A, b, block_size, 'cellMajorOrder', true, ....
                'tolerance', tol, 'block_size', 2, 'use_drs', false);
            test.assertEqual(x, ref, 'AbsTol', 10*tol)
            test.assertFalse(err > tol);
        end
        
        function CPRBlockTest(test)
            [A, b, ref] = test.getBlockMatrix();
            tol = 1e-8;
            block_size = 2;
            [x, err] = callAMGCL_cpr(A, b, block_size, 'cellMajorOrder', true, ....
                'tolerance', tol, 'block_size', 2, 'use_drs', false);
            test.assertEqual(x, ref, 'AbsTol', 10*tol)
            test.assertFalse(err > tol);
        end
        
        function CPRDRSScalarTest(test)
            [A, b, ref] = test.getBlockMatrix(1);
            tol = 1e-8;
            block_size = 2;
            [x, err] = callAMGCL_cpr(A, b, block_size, 'cellMajorOrder', true, ....
                'tolerance', tol, 'block_size', 2, 'use_drs', true);
            test.assertEqual(x, ref, 'AbsTol', 10*tol)
            test.assertFalse(err > tol);
        end
        
        function CPRDRSBlockTest(test)
            [A, b, ref] = test.getBlockMatrix();
            tol = 1e-8;
            block_size = 2;
            [x, err] = callAMGCL_cpr(A, b, block_size, 'cellMajorOrder', true, ....
                'tolerance', tol, 'block_size', 2, 'use_drs', true);
            test.assertEqual(x, ref, 'AbsTol', 10*tol)
            test.assertFalse(err > tol);
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

