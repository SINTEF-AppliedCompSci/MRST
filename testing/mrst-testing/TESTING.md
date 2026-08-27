# MRST Pre-Release Testing Workflow

This document describes the three-tier pre-release testing system introduced
in MRST and explains how module authors can integrate their code with it.

---

## Overview

| Tier | Name                     | When to run           | Tooling                              |
|------|--------------------------|-----------------------|--------------------------------------|
| 1    | Unit tests               | Every push / PR       | `matlab.unittest.TestCase`           |
| 2    | Integration / regression | Local / ad hoc runs   | `matlab.unittest` + `MRSTRegressionTest` |
| 3    | Example smoke tests      | Local / ad hoc runs   | `MRSTExampleTests` (parameterised)   |

All three tiers are orchestrated by a single entry point:

```matlab
runPreReleaseTests('tier', [1 2 3])
```

---

## Running Tests Locally

```matlab
% Fast pre-commit check (Tier 1 only)
startup();
mrstModule add ad-unittest;
runPreReleaseTests('tier', 1);

% Nightly / pre-merge check (Tier 1 + 2)
runPreReleaseTests('tier', [1 2]);

% Full pre-release for specific modules
runPreReleaseTests('tier', [1 2 3], 'modules', {'ad-blackoil', 'co2lab'});

% Run slow examples that CI skips by default
runPreReleaseTests('tier', 3, 'excludeTags', {});

% Audit module coverage before a release
mrstModuleTestStatus();
```

JUnit XML output can be requested with `'writeXML', true`.

---

## GitHub Actions CI

The current CI workflow only runs Tier 1 unit tests:

| File                                    | Trigger                            | Tiers |
|-----------------------------------------|------------------------------------|-------|
| `.github/workflows/unit-tests.yml`      | Every push / PR to `main`, `release-*` | 1 |

Tier 2 and Tier 3 are still available through `runPreReleaseTests`, but they are
not scheduled by default CI in the current branch.

---

## Module Author Interface

Any MRST module *may* provide the following artefacts to participate in the
testing system.  None are mandatory — the system degrades gracefully when they
are absent.

### 1. Unit test classes (`tests/unitTests/Test*.m` or `tests/unitTests/*Test.m`)

Place `matlab.unittest.TestCase` subclasses inside a `tests/unitTests/`
directory at the root of your module. Files matching `Test*.m` or `*Test.m`
are auto-discovered by `getUnitTestSuiteMRST`.

Copy the template from:
```
autodiff/ad-unittest/setup/moduleTestTemplate/TestMyModule.m
```

Example:

```matlab
classdef TestMyFluid < matlab.unittest.TestCase
    methods (TestClassSetup)
        function loadModules(tc)
            mrstModule add ad-core my-module
        end
    end
    methods (Test)
        function densityAtSurface(tc)
            fluid = initSimpleADIFluid();
            tc.verifyGreaterThan(fluid.rhoOS, 0);
        end
    end
end
```

### 2. Regression test registration (`getModuleRegressionTests.m`)

Place a function `getModuleRegressionTests.m` at the root of your module.
It must return a cell array of `MRSTRegressionTest` objects. These are run as
part of Tier 2.

Copy the template from:
```
testing/mrst-testing/setup/moduleTestTemplate/getModuleRegressionTests.m
```

Each regression test should be easy to read and run:

- define it with `MRSTRegressionTest`
- keep the setup in a plain function
- use `packSimulationProblem`/`copyPackedProblem` for packed-problem setups
- let the shared comparison helper handle states, well solutions, and reports
- if you need a setup struct, use `processMRSTRegressionSetupInput` and
  `packRegressionSetup`

### 3. Example skip list (`getSkippedExamples.m`)

Place a function `getSkippedExamples.m` at the root of your module to
declare examples that should never run in automated testing (e.g. examples
that depend on proprietary datasets you cannot distribute).

Prefer per-file annotations (see below) for single-example exclusions.

### 4. Per-file `MRST_TEST_OPTIONS` annotations

Add a structured comment block near the top of any example file
(within the first 50 lines):

```matlab
%{
MRST_TEST_OPTIONS
timeout:     120        % seconds; 0 = no limit (default)
interactive: false      % true = skip in all automated runs
tags:        slow       % comma-separated: slow, data-required, gpu, ...
%}
```

**Available tags:**

| Tag             | Meaning                                                         |
|-----------------|-----------------------------------------------------------------|
| `slow`          | Example takes several minutes; excluded from default CI run    |
| `data-required` | Requires external datasets not bundled with MRST               |
| `interactive`   | Requires GUI / user interaction                                |
| `gpu`           | Requires a GPU                                                 |

The template file `myExampleTemplate.m` shows complete usage:
```
autodiff/ad-unittest/setup/moduleTestTemplate/myExampleTemplate.m
```

---

## Regenerating Regression Baselines

After an intentional API change that alters simulation output, baselines must
be regenerated.  Use the `RegressionTest` workflow:

```matlab
% Delete existing baseline and create a new one
rt = RegressionTest('my_test_case', 'group', 'mymodule');
rt.runRegressionTest('deletePrevious', true);
```

Commit the updated reference files from `ad-unittest/output/` to the
repository.

---

## Auditing Test Coverage

```matlab
mrstModuleTestStatus()
```

Prints a table showing, for each module:
- Number of unit test files discovered
- Whether a `getModuleRegressionTests.m` is present
- Whether a `getSkippedExamples.m` is present
- Total examples and how many carry `MRST_TEST_OPTIONS` annotations

Use this before a release to identify modules with no testing coverage.
