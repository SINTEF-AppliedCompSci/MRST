classdef CoreADAndUtilitiesIntegrationTest < matlab.unittest.TestCase
    methods
        function test = CoreADAndUtilitiesIntegrationTest()
            mrstModule reset
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
