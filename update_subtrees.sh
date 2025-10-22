#!/bin/bash

# This repository is a repo that includes all released MRST parts using git subtree.
# Running this script pulls all subtrees.
# Intentionally kept simple. Note that main/master depends on when repo was created.
#
# To add new subtree:
# git subtree add --prefix {local path} {repo URL} {remote branch}
# To pull:
# git subtree pull -m "Update subtree" --prefix {local_path} {remote} {remote branch}

# Single module repos (third party)
git subtree pull -m "Update subtree" --prefix modules/ad-micp ad-micp main
git subtree pull -m "Update subtree" --prefix modules/dual-continuum-mech dual-continuum-mech master
git subtree pull -m "Update subtree" --prefix modules/dual-porosity dual-porosity master
git subtree pull -m "Update subtree" --prefix modules/dual-porosity-permeability dual-porosity-permeability main
git subtree pull -m "Update subtree" --prefix modules/fv-unsat fv-unsat master
git subtree pull -m "Update subtree" --prefix modules/fvbiot fvbiot master
git subtree pull -m "Update subtree" --prefix modules/hwu-fractures hwu-fractures master
git subtree pull -m "Update subtree" --prefix modules/nwm nwm master
git subtree pull -m "Update subtree" --prefix modules/re-mpfa re-mpfa master
git subtree pull -m "Update subtree" --prefix modules/shale shale main
git subtree pull -m "Update subtree" --prefix modules/upr upr master
