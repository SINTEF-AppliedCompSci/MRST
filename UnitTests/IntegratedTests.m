numTests = 3;
isPassed = zeros(numTests,1);
isPassed(1) = testEDFMNFR();
isPassed(2) = testpEDFMNFR();
isPassed(3) = testDiffusion();

fprintf('%d Tests passed out of %d \n',sum(isPassed),numTests);