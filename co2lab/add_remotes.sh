#!/bin/bash
# Commands to add all remotes
# Intentionally kept simple.

# Multi-module repos and core
git remote add core ssh://git@bitbucket.org/mrst/mrst-core.git
git remote add autodiff ssh://git@bitbucket.org/mrst/mrst-autodiff.git
git remote add model-io ssh://git@bitbucket.org/mrst/mrst-model-io.git
git remote add multiscale ssh://git@bitbucket.org/mrst/mrst-multiscale.git
git remote add visualization ssh://git@bitbucket.org/mrst/mrst-visualization.git
git remote add solvers ssh://git@bitbucket.org/mrst/mrst-solvers.git

# Single module repos
git remote add co2lab ssh://git@bitbucket.org/mrst/mrst-co2lab.git
git remote add book ssh://git@bitbucket.org/mrst/mrst-book.git
git remote add dg ssh://git@bitbucket.org/mrst/mrst-dg.git
git remote add nfvm ssh://git@bitbucket.org/mrst/nfvm.git
git remote add mpsaw ssh://git@bitbucket.org/mrst/mpsaw.git
git remote add ensemble ssh://git@bitbucket.org/mrst/mrst-ensemble.git
git remote add network-models ssh://git@bitbucket.org/mrst/network-models.git
git remote add geothermal ssh://git@bitbucket.org/mrst/geothermal.git
git remote add static-modeling ssh://git@bitbucket.org/mrst/static-modeling.git

# Single module repos (third party)
git remote add ad-micp ssh://git@github.com/daavid00/ad-micp.git
git remote add domain-decomposition ssh://git@bitbucket.org/mrst/mrst-domain-decomposition.git
git remote add dual-continuum-mech ssh://git@bitbucket.org/mashworth92/dual-continuum-mech.git
git remote add dual-porosity ssh://git@bitbucket.org/rafael_march_carbonates_hw/mrst-dual-porosity.git
git remote add dual-porosity-permeability ssh://git@github.com/nikolai-andrianov/DPDP-MRST.git
git remote add fv-unsat ssh://git@github.com/jhabriel/fv-unsat.git
git remote add fvbiot ssh://git@github.com/pmgbergen/fvbiot.git
git remote add geochemistry ssh://git@bitbucket.org/mrst/matlab-geochemistry.git
git remote add hwu-fractures ssh://git@bitbucket.org/HWUCarbonates/mrst-hwu-fractures.git
git remote add nwm ssh://git@bitbucket.org/LinZhao9/nwm.git
git remote add re-mpfa ssh://git@github.com/jhabriel/RE-MPFA.git
git remote add shale ssh://git@github.com/UnconvRS/shale.git
git remote add upr ssh://git@github.com/rbe051/UPR.git
