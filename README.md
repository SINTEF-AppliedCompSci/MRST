# MRST: The Matlab Reservoir Simulation Toolbox

## What is this?
This is the official GitHub mirror of the Matlab Reservoir Simulation Toolbox (MRST), an open source toolbox for simulation of flow, mechanics and transport in porous media developed at [group for Computational Geosciences](https://www.sintef.no/en/digital/applied-mathematics/computational-geoscience/) at [SINTEF Digital](https://www.sintef.no/digital/). We only recommend getting MRST from the repositories for users who are interested in following the development actively and tolerate the occasional bug or broken feature. For most users, the [MRST releases](https://www.sintef.no/projectweb/mrst/download/) released twice each year are a better option, which comes pre-packaged and tested.

For more details on MRST, please see our website at [www.mrst.no](http://www.mrst.no) and the [MRST repositories on Bitbucket](https://www.bitbucket.org/mrst).

## Getting started with MRST
Clone this repository:

```bash
git clone https://github.com/SINTEF-COMG/MRST.git
```

Once the repository has been cloned, navigate MATLAB or GNU Octave to the checked out repository

```matlab
cd MRST;
run startup.m
```

This should initialize MRST and produce a welcome message:
```
Welcome to the MATLAB Reservoir Simulation Toolbox (MRST)!
You are using the release version 2022b. To download other versions of MRST
and view examples and relevant publications, please visit www.mrst.no

Useful commands for getting started:
 - List all introductory examples:   mrstExamples()
 - List all modules:                 mrstPath('list')
 - Load modules using GUI:           mrstModule('gui')
 - Explore all available data sets:  mrstDatasetGUI()
 - List examples of a module:        mrstExamples('ad-blackoil')
 - Explore modules and publications: mrstExploreModules()
 - Show all examples in all modules: mrstExamples('all')
 - Show settings for MRST:           mrstSettings()
 - Display this message:             mrstStartupMessage()

For assistance and discussions about MRST, please visit our mailing list at
   www.sintef.no/projectweb/mrst/forum/ (sintef-mrst@googlegroups.com)
For some common queries, see our FAQ: www.sintef.no/projectweb/mrst/faq/
```

You can verify your installation by running a basic tutorial that produces a plot:
```matlab
flowSolverTutorial1
```

Note that MRST takes care not to modify files outside the checked out directory. For this reason, you will have to re-run the `startup.m` function if you restart MATLAB/Octave.

## Additional resources

To learn more about MRST, we suggest that you check out the [documentation that includes two free open-access books](https://www.sintef.no/projectweb/mrst/documentation/).

## Access to repositories
MRST is currently organized into several repositories. These are added to this repository using `git subtree` for your convenience together with a `startup.m` script that can be run with MATLAB and GNU Octave to load these files into your path. This is primarily intended as a mirror for getting access to the development version through a single git clone. If you would rather work directly with individual repositories, please follow the instructions at the [MRST core wiki at Bitbucket](https://bitbucket.org/mrst/mrst-core/wiki/Home).


| Name                  | Description                                                                                                                                  | URL                                                   |
|--------------------   |-------------------------------------------------------------------------------------------------------------------------------------------   |----------------------------------------------------   |
| core                  |  Core functionality for MRST                                                                                                                 | https://bitbucket.org/mrst/mrst-core                  |
| autodiff              |  Solvers based on automatic differentiation: Black-oil, compressible flow, compositional and the AD framework                                | https://bitbucket.org/mrst/mrst-autodiff              |
| solvers               | Different solvers without automatic differentiation: Consistent discretizations, fractured media and geomechanics.                           | https://bitbucket.org/mrst/mrst-solvers               |
| co2lab                | Numerical CO2 laboratory: Solvers and utilities for large-scale CO2 injection.                                                               | https://bitbucket.org/mrst/mrst-co2lab                |
| model-io              | Different utilities and modules for input and output to MRST.                                                                                | https://bitbucket.org/mrst/mrst-model-io              |
| multiscale            | Solvers and routines for scale transitions, including multiscale and upscaling implementations.                                              | https://bitbucket.org/mrst/mrst-multiscale            |
| visualization         | Modules with a focus on visualization and graphical user interfaces.                                                                         | https://bitbucket.org/mrst/mrst-visualization         |
| thirdparty-modules    | Third-party modules for MRST. This repository contains a set of git submodules which link to other git repositories with useful modules.     | https://bitbucket.org/mrst/mrst-thirdparty-modules    |

The remaining modules are individual repositories and can be found in the `add_remotes.sh` file.

## Updating the subtree

Normally you can pull this repository to get updates to MRST. There is some delay between pushes to the individual module repositories and this repository as it is a manual process. 
### Advanced usage - updating subtrees manually
If you want to update manually, or if you have a local fork of this repository and you want to update the files, you should first do the following:
```bash
./add_remotes.sh
```
This needs only to be done once for a given clone of MRST. Once the remotes are set up and the files are executable you can update the subtrees :

```
./update_subtrees.sh
```
Note that updating the subtrees create merge commits. If you want to go back in sync with the repository (e.g. to do a git pull) you can run the following. 

**This will delete any local changes you have and reset it like a fresh clone had been performed.**
``` bash
git fetch origin
git reset --hard origin/master
```