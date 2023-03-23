numTests = 6;
isPassed = zeros(numTests,1);
isPassed(1) = testPEDFM_EDFM_explicit();
isPassed(2) = testEDFMNFR();
isPassed(3) = testpEDFMNFR();
isPassed(4) = testDiffusion();
isPassed(5) = testSorption();
isPassed(6) = testGangi();

fprintf('%d Tests passed out of %d \n',sum(isPassed),numTests);