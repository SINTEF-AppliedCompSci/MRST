![MRST Logo](mrst.png)
## What is this?
This is the official GitHub repository of the Matlab Reservoir Simulation Toolbox (MRST), an open source toolbox for simulation of flow, mechanics and transport in porous media developed at [group for Applied Computational Science](https://www.sintef.no/en/digital/departments-new/department-of-mathematics-and-cybernetics/research-group-applied-computational-science/) at [SINTEF Digital](https://www.sintef.no/digital/). We only recommend getting MRST from the repositories for users who are interested in following the development actively and tolerate the occasional bug or broken feature. For most users, the [MRST releases](https://www.sintef.no/projectweb/mrst/download/) released twice each year are a better option, which comes pre-packaged and tested.

For more details on MRST, please see our website at [www.mrst.no](http://www.mrst.no). Starting with MRST 2024b, development of MRST occurs in this repository. The following GitHub features are enabled if you want to get involved:

 - [**Discussions**](https://github.com/SINTEF-AppliedCompSci/MRST/discussions) can be used for questions/discussions.
 - [**Issues**](https://github.com/SINTEF-AppliedCompSci/MRST/issues) for filing bugs or feature requests.
 - [**Pull requests**](https://github.com/SINTEF-AppliedCompSci/MRST/pulls) for feature contributions.

## Getting started with git version of MRST
Clone this repository:

```bash
git clone git@github.com:SINTEF-AppliedCompSci/MRST.git
```

Once the repository has been cloned, navigate your install of MATLAB or GNU Octave to the checked out repository

```matlab
cd MRST;
run startup.m
```

This should initialize MRST and produce a welcome message:
```
Welcome to the MATLAB Reservoir Simulation Toolbox (MRST)!
You are using the release version 2024b. To download other versions of MRST
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

## Updating subtrees for third party modules

Some of the external modules developed for MRST are individual repositories and can be found in the `add_remotes.sh` file. We update these from time to time and the commits get added to the main repository.

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
