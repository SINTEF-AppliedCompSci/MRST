# Testing MRST

MRST's automated tests are organized into three tiers, all run through a
single entry point, `runMRSTTests`.

| Tier         | What it covers                                             | When it runs                        |
|--------------|-------------------------------------------------------------|--------------------------------------|
| `unit`       | Fast, isolated tests: `core/` + autodiff (`ad-unittest`)     | Every push/PR (`.github/workflows/unit-tests.yml`) |
| `regression` | Heavier scenario/simulation tests (autodiff's `test_sim`)    | Manual only (`.github/workflows/regression-tests.yml`) |
| `examples`   | Runs every example script in a module end-to-end            | Manual only, per module (`.github/workflows/example-smoke-tests.yml`) |

`autodiff/test-suite` is a separate, older regression/benchmark harness and
is **not** part of this system.

## Running tests locally

After `startup`, from the MATLAB command line:

```matlab
% Unit tests: core/ + autodiff
runMRSTTests('unit');

% Regression tests: autodiff scenario suite
runMRSTTests('regression');

% Example smoke tests: one module
runMRSTTests('examples', 'module', 'ad-blackoil');

% Example smoke tests: every registered module (slow -- ad hoc use only)
runMRSTTests('examples', 'module', 'all');
```

Useful options (see `help runMRSTTests` for the full list):
- `'writeXML', true` -- write JUnit XML reports (what CI uses).
- `'errorOnFailure', true` -- raise a MATLAB error if any test failed,
  instead of just printing the summary (what CI uses, for a non-zero exit
  code).

Every run ends with a command-line summary: a per-suite pass/fail/incomplete
count, an explicit list naming every failed test, and an overall total.
Start there when a run fails -- it tells you exactly which suite and test
to look at without scrolling back through the full run log.

## Adding a new core unit test

`core/` has no automatic test discovery, because some files under `core/`
that look like tests by name are actually old demo/plotting scripts, not
real `matlab.unittest` tests, and blind folder scanning would trip over
those. Instead:

1. Write a `matlab.unittest`-conformant test file near the code it covers
   -- either function-based:
   ```matlab
   function tests = myThingTest
       tests = functiontests(localfunctions);
   end

   function testSomeBehavior(t)
       verifyEqual(t, myThing(2), 4);
   end
   ```
   or classdef-based (`classdef myThingTest < matlab.unittest.TestCase`).
   See `core/utils/gridtools/tests/gridToolsTest.m` or
   `core/utils/equil/simpleEquilibriumTest.m` for worked examples.
2. Add its full path to the explicit list in
   `core/utils/testing/getCoreUnitTestSuiteMRST.m`.

## Adding a new autodiff unit or regression test

Unlike `core/`, these folders under `autodiff/ad-unittest/` *are*
auto-discovered (`matlab.unittest.TestSuite.fromFolder`) -- no registration
step needed, just place a conformant test file in the right folder:

- `autodiff/ad-unittest/test_models/` or `test_utils/` -- unit tier.
- `autodiff/ad-unittest/test_sim/` -- regression tier.

## Adding a new example to a module

This is deliberately zero-effort: drop your `.m` script under
`<module>/examples/` (any subfolder is fine except `utils` or `helpers`,
which are excluded). `mrstExamples` discovers it automatically, and it
becomes a smoke-test case the next time `runMRSTTests('examples', 'module', '<module>')`
runs for that module -- no code change, no list to edit.

The one thing you may need to do: if the new example is slow, opens a GUI,
or launches external MATLAB sessions (like the existing skipped examples
`runNorneExample`, `preprocessDiagnosticsEgg`, `ensembleGUIForEgg`), add its
filename to the single skip-list in `getSkippedTests()` inside
`autodiff/ad-unittest/MRSTExampleTests.m`. Don't create a second skip-list
elsewhere -- this is the one place the example-smoke tier checks.

## CI workflows

- `.github/workflows/unit-tests.yml` -- runs automatically on push/PR.
- `.github/workflows/regression-tests.yml` -- manual only: open the
  Actions tab, select "Regression Tests", and click "Run workflow".
- `.github/workflows/example-smoke-tests.yml` -- manual only, same as
  above, with a "module" input (a module name, or `all` for every module).

All three run on real MATLAB (via MathWorks' `matlab-actions`, free for
public repositories), not Octave -- the existing test suites are built on
`matlab.unittest`, which Octave does not implement.
