# The effect of Alzheimer's disease and its progression on pyramidal cell gain and connectivity underlying mismatch negativity responses
_This code accompanies Lanskey et al. (2024). The code runs on Matlab 2020b and uses SPM12 v7771, which is freely available:
* [Statistic Parametric Mapping](http://www.fil.ion.ucl.ac.uk/spm/)

## Analysis steps
The analysis is run in the following order:
1. Preprocessing, using the folder '1_preprocessing'
2. Sensor space analysis, folder '2_sensor space'
3. First-level DCM inversion, folder '3_dcm'
4. Second-level Parameteric Empirical Bayes, folder 4_peb

### Preprocessing ('/1_preprocessing')
Preprocessing is completed using the following scripts:

```
ntad_mmn_preprocessing_v2_jlrep10
```
This routine uses a standard preprocessing pipeline developed for the NTAD study. Preprocessing is completed for all participants and all sessions (baseline, two week, follow up)

```
sensor_combinebc_v1
```
This script combines planar gradiometers for the sensor-level analyses

```
coregister_v1
```
This script coregisters each subjects T1-weighted MRI to their MEG files for DCM source inversion

### Sensor space analysis ('/2_sensor_space')
Sensor space analyses uses the following scripts

```
sensor_BF_v1
```
This script runs the sensor space t-tests and ANCOVAs for the baseline (controls vs AD/MCI) and longitudinal (AD/MCI baseline vs follow up) comparisons and generates mismatch negativity waveform plots.


```
sensor_bltw_v1
```
This script calculates the absolute intraclass correlation of the mismatch negativity amplitude for the baseline and two-week MEG scans

### First-level DCM inversion ('/3_dcm')
DCMs are inverted for baseline and follow-up scans using the following:

```
dcm_cmc_BL_v1
```
This routine inverts the full DCM for each individual at baseline and checks model fits, then runs a bayesian model reduction and compares models to identify how the parietal node is connected to other regions


```
dcm_cmc_BL_v1
```
This routine inverts the full DCM for each individual at follow up and checks model fits, then runs a bayesian model reduction and compares models to identify how the parietal node is connected

### Second-level PEB ('/4_peb)
This script compares parameters of interest at baseline (controls versus AD/MCI) and longitudinally (AD/MCI baseline versus AD/MCI follow up)

## Custom functions
* `cmc_environment` - This function sets the environment and folder structures
* `mmn_amp` - calculates the MMN amplitude
* `dcm_cmc_gen_model_space_v1` - generates alternate model structures (with different parietal connectivity in A and C matrices)
* `dcm_cmc_gen_v1` - specifies DCM options
* `dcm_cmc_gen_expo_v1` - specifies DCM options with an exponential repetition effect only (i.e. no phasic repition effect)
* `plot_bmc_F_jl` - Generates a plot with differential free energies across models (top panel) and posterior model probabilities (lower panel)
* `spm_dcm_peb_review_fig_jl` - Adapted from the SPM script 'spm_dcm_peb_review' to generate PEB figures
