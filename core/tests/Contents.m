% TESTS
%   MATLAB unittest coverage for core MRST primitives and workflows.
%   Tests are organized into two tiers for use with runPreReleaseTests.m:
%     - Tier 1 (Unit Tests): Fast tests in unitTests/ subdirectory.
%     - Tier 2 (Integration Tests): Slower tests in integrationTests/ subdirectory.
%
% Subdirectories
%   unitTests                   - Tier 1 unit tests (fast, < 1 min each)
%   integrationTests            - Tier 2 integration tests (< 10 min each)
%
% Files (see subdirectories for test classes)
%   runTestsCore                 - [DEPRECATED] Use runPreReleaseTests.m instead.
%
% Test Classes in unitTests/
%   CoreGridGeometryTest         - Grid, topology, and geometry tests (tier1)
%   CoreRockFlowPrimitivesTest   - Rock/state/initialization tests (tier1)
%   CoreADAndUtilitiesTest       - AD/utilities tests (tier1)
%
% Test Classes in integrationTests/
%   CoreCornerPointProcessingTest     - Corner-point processing tests (tier2)
%   CoreRockFlowPrimitivesIntegrationTest - Flow setup tests (tier2)
%   CoreADAndUtilitiesIntegrationTest     - Module plumbing tests (tier2)
%
% NOTE:
%   The folder structure follows the runPreReleaseTests.m discovery conventions:
%     - tests/unitTests/ for Tier 1 unit tests
%     - tests/integrationTests/ for Tier 2 integration tests
%   To run all tests: runPreReleaseTests('tier', [1 2], 'modules', {'ad-core'})
%   To run only tier 1: runPreReleaseTests('tier', 1, 'modules', {'ad-core'})

