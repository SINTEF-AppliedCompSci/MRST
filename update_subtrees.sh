#!/bin/bash

# This repository is a repo that includes all released MRST parts using git subtree.
# Running this script pulls all subtrees.
# Intentionally kept simple. Note that main/master depends on when repo was created.
#
# To add new subtree:
# git subtree add --prefix {local path} {repo URL} {remote branch}
# To pull:
# git subtree pull -m "Update subtree" --prefix {local_path} {remote} {remote branch}

# Multi-module repos and core
git subtree pull -m "Update subtree" --prefix core core master
git subtree pull -m "Update subtree" --prefix autodiff autodiff master
git subtree pull -m "Update subtree" --prefix model-io model-io master
git subtree pull -m "Update subtree" --prefix multiscale multiscale master
git subtree pull -m "Update subtree" --prefix visualization visualization master
git subtree pull -m "Update subtree" --prefix solvers solvers master

# co2lab is special - contains both modules and a deprecated folder
git subtree pull -m "Update subtree" --prefix co2lab co2lab master

# Single module repos
git subtree pull -m "Update subtree" --prefix modules/network-models network-models master
git subtree pull -m "Update subtree" --prefix modules/book book master
git subtree pull -m "Update subtree" --prefix modules/dg dg master
git subtree pull -m "Update subtree" --prefix modules/nfvm nfvm master
git subtree pull -m "Update subtree" --prefix modules/mpsaw mpsaw master
git subtree pull -m "Update subtree" --prefix modules/ensemble ensemble master
git subtree pull -m "Update subtree" --prefix modules/geothermal geothermal master
git subtree pull -m "Update subtree" --prefix modules/static-modeling static-modeling master

# Single module repos (third party)
git subtree pull -m "Update subtree" --prefix modules/ad-micp ad-micp main
git subtree pull -m "Update subtree" --prefix modules/domain-decomposition domain-decomposition master
git subtree pull -m "Update subtree" --prefix modules/dual-continuum-mech dual-continuum-mech master
git subtree pull -m "Update subtree" --prefix modules/dual-porosity dual-porosity master
git subtree pull -m "Update subtree" --prefix modules/dual-porosity-permeability dual-porosity-permeability main
git subtree pull -m "Update subtree" --prefix modules/fv-unsat fv-unsat master
git subtree pull -m "Update subtree" --prefix modules/fvbiot fvbiot master
git subtree pull -m "Update subtree" --prefix modules/geochemistry geochemistry master
git subtree pull -m "Update subtree" --prefix modules/hwu-fractures hwu-fractures master
git subtree pull -m "Update subtree" --prefix modules/nwm nwm master
git subtree pull -m "Update subtree" --prefix modules/re-mpfa re-mpfa master
git subtree pull -m "Update subtree" --prefix modules/shale shale main
git subtree pull -m "Update subtree" --prefix modules/upr upr master
