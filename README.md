[![Build Status](https://github.com/daavid00/ad-micp/actions/workflows/CI.yml/badge.svg)](https://github.com/daavid00/ad-micp/actions/workflows/CI.yml)[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# ad-micp: A module to study CO2 leakage remediation by microbially induced calcite precipitation (MICP)

<img src="micp-gif.gif" width="830" height="400">

## Description
MICP is a new and sustainable technology which utilizes biochemical
processes to create barriers by calcium carbonate cementation; therefore,
this technology has a potential to be used for sealing leakage zones in
geological formations. We have implemented a mathematical model for MICP
suitable for field-scale studies. Further information on the mathematical
model can be found in the two publications [A, B] in the "Cite" part of
this README.

The numerical examples accompanying this module are:
* runMICP1DCase.m and
* runMICP3DCase.m.

In addition, this module includes scripts to run the studies in [A, B], see the [_publications folder_](https://github.com/daavid00/ad-micp/tree/main/publications).

The ad-micp module is compatible with the (freely available) MATLAB
Reservoir Simulation Toolbox (MRST) provided by SINTEF Digital, see
http://www.sintef.no/projectweb/mrst/.

See the [_Wiki_](https://github.com/daavid00/ad-micp/wiki) for an extended description, information regarding changes
respect to previous releases, and using this module in macOS with GNU Octave.

## Requirements
* MRST (Tested version: mrst-2024b)
* MATLAB (Tested version: R2023a) or GNU Octave (Tested version: 9.2.0)

## Installation (MATLAB/GNU Octave)
The ad-micp module is included in MRST, click [_this link_](https://github.com/SINTEF-AppliedCompSci/MRST/tree/main/modules/ad-micp) 
to see the version of ad-micp included.

If you are interested in using the ad-micp module from the edge version, this can be achieved by the following steps:

```bash
# Clone MRST
git clone https://github.com/SINTEF-AppliedCompSci/MRST.git
# Get inside the modules folder
cd MRST/modules
# Remove the ad-micp folder from the MRST folder
rm -rf ad-micp
# Clone the ad-micp repo
git clone https://github.com/daavid00/ad-micp.git
``` 

Then you can try to run the tests, e.g., the [_test_runMICP1DCase.m_](https://github.com/daavid00/ad-micp/blob/main/tests/test_runMICP1DCase.m) using octave (octave can be installed by executing `brew install octave`):

```bash
octave tests/test_runMICP1DCase.m
```

## Getting started
See https://www.youtube.com/watch?v=nvz3bV4QgxM.

## Cite
If you use ad-micp to write a scientific publication, please cite one of
the following papers:
* [A] Landa-Marbán, D., Tveit, S., Kumar, K., Gasda, S.E., 2021.
Practical approaches to study microbially induced calcite precipitation
at the field scale. Int. J. Greenh. Gas Control 106, 103256.
https://doi.org/10.1016/j.ijggc.2021.103256.
* [B] Landa-Marbán, D., Kumar, K., Tveit, S., Gasda, S.E.
Numerical studies of CO2 leakage remediation by micp-based plugging
technology. In: Røkke, N.A. and Knuutila, H.K. (Eds) Short Papers from the
11th International Trondheim CCS conference, ISBN: 978-82-536-1714-5,
284-290.

## Contact
David Landa-Marbán (dmar@norceresearch.no).
