# The Soil Moisture Index - SMI program v2.0.5

This repository contains the soil moisture index (SMI) Fortran program developed at the Dept. Computational Hydrosystems at the Helmholtz Centre for Environmental Research - UFZ.

The SMI comes with a [LICENSE][1] agreement, this includes also the GNU Lesser General Public License.

**Please note**: The GitLab repository grants read access to the code.
If you like to contribute to the code, please contact stephan.thober@ufz.de.

## Installation

Installation instructions can be found in [INSTALL][2** for Windows, MacOS, and GNU/Linux distributions.

## Release Notes

**Version 2.0.5**

- Removed dependency on numerical recipes

**Version 2.0.4**

- includes flexible latlon coordinates (1D and 2D)

**Version 2.0.3**

- BUG FIX: openmp loop for SMI inversion now works correctly

**Version 2.0.2**

- BUG FIX: Check for monthly data covering multiple years is now correct.

**Version 2.0.1**

- implemented openmp parallelization for estimation of SMI values

**Version 2.0.**

- consistent handling of leap days, no kde estimated for leap days,
  instead data of March 1st is used.
- soil moisture and kernel widths for kde are stored in one netcdf
  file called cdf_info.nc.
- added period object in mo_global_variables.
- mask file is now optional

## Bug fixes

- kde is correctly applied for sub-annual data

## Usage

The SMI code uses a kernel density estimator with a gaussian kernel for transforming soil
moisture values to the soil moisture index (SMI). It has several
features that are described in the namelist file main.dat.

## Cite as

Please refer to the main model by citing Samaniego et al. (2013) and Samaniego et al. (2018):

- Samaniego et al. (2013), "Implications of Parameter Uncertainty on Soil Moisture Drought Analysis in Germany", J Hydrometeor, 2013 vol. 14 (1) pp. 47-68. http://journals.ametsoc.org/doi/abs/10.1175/JHM-D-12-075.1

- Samaniego et al. (2018), "Anthropogenic warming exacerbates European soil moisture droughts", Nature Climate change, 2018 pp. 1-9. http://dx.doi.org/10.1038/s41558-018-0138-5

[1]: LICENSE
[2]: INSTALL.md
