SnPM13: Bugs & Fixes

This file describes the bugs that have been reported, along with the appropriate fixes. 
Updated versions of appropriate SnPM functions are available from in the snpm13_updates: 
http://warwick.ac.uk/snpm/distribution/snpm13_updates

--- SnPM 13.1.06 ---
 * fix: ANOVA between group "all zero"/"all equal" options were flipped in the menu of the matlabbatch.
 * fix: snpm_cp can now handle XYZ locations in a 1D space.
 * Version 13.1.06 (`snpm`)

--- SnPM 13.1.05 ---
 * fix: two-sample t-test with nscans>12 was errored due to call to undefined variable `nPerm` (bug introduced in previous release: 13.1.4) (`snpm_pi_TwoSampT`)
 * fix: warning on size of `SnPM_ST.mat` was never raised
 * Add seed shuffling for paired two-sample test modules (`snpm_pi_PairT`, `snpm_pi_PairTrand`)
 * Version 13.1.05 (`snpm`)

--- SnPM 13.1.04 ---
* snpm
SnPM version 13.1.04

* snpm_bch_ui, snpm_cp, snpm_ui
Fix: implicit masking for non-float datatypes

* snpm_pi_ANOVAbtwnS, snpm_pi_ANOVAwithinS, snpm_pi_Corr, snpm_pi_Corr1S, 
snpm_pi_OneSampT, snpm_pi_TwoSampPairT, snpm_pi_TwoSampT, snpm_pi_TwoSampTss
Update syntax to set the seed of the random number generator (rng) 

* snpm_defaults
New SnPM global parameter shuffle_seed. If true, the rng uses a 'shuffle' seed providing different results for each run. If false Matlab's 'default' seed or any seed defined by the user will be used .

* snpm_pi_TwoSampT
Allows MC permutation (with reptitions) for large enough nPerms 

* snpm_pp, snpm_ui
Do not display figures in command-line mode

* snpm, snpm_check_nperm, snpm_combo, snpm_cp, snpm_init, snpm_pi_ANOVAbtwnS,
snpm_pi_ANOVAwithinS, snpm_pi_Corr, snpm_pi_Corr1S, snpm_pi_OneSampT, 
snpm_pi_PairT, snpm_pi_PairTrand, snpm_pi_TwoSampPairT, snpm_pi_TwoSampT, 
snpm_pi_TwoSampTss, snpm_pp, snpm_t2z, snpm_uc_FDR, snpm_ui, spm_append_96
Add identifiers to errors and warnings

* snpm_cp
Fix: Do not assume statistic is T when computing cluster forming threshold from P-value

* generic_test_snpm, skeleton_class, test_ANOVAbetween, 
test_multisub_withinsubanova, test_multisubpaired2cond, 
test_multisubsimpleregression, test_oneSample, test_onesub_regression, 
test_onesub_twocondrepl, test_twoSample, test_twosample_twocond
Test in command-line mode with no shuffling for the random seed, using new rng
 syntax.

--- Updates from SnPM 13.1.03 ---
* snpm
SnPM version 13.1.03

* snpm_pi_PairT
Now allows more than 52 subjects, as previously that would generate a "Maximum variable size allowed by the program" error message.

--- Updates from SnPM 13.1.02 ---
* snpm
SnPM version 13.1.02

* snpm_pp, snpm_combo_pp
Update contrast display for compatibility with Matlab R2014b.

* generic_test_snpm
Updates affecting tests only.

--- Updates from SnPM 13.1.01 ---
* snpm
SnPM version 13.1.01

* snpm_pp
Check that 'Locs_vox' is integer according to a pre-defined tolerance level (as in snpm_cp).

--- Updates from SnPM13.1.00 ---
* snpm_cp, snpm_pp, generic_test_snpm, test_oneSample
Cluster-forming threshold defined by a P-value are now extended to pseudo-T statistic (with variance smoothing). The Z-score corresponding to the selected P-value is used as an approximation of pseudo-T values.

--- Updates from SnPM13.0.14 ---
* snpm_bch_ui_Corr, snpm_bch_ui_Corr1S, snpm_pi_Corr, snpm_pi_Corr1S, test_multisubsimpleregression, test_onesub_regression
Introduce nuisance covariates in simple-subject and multi-subject regressions.

--- Updates from SnPM13.0.13 ---
* snpm_combo_pp
Fix: Set bNeg to 0 for F-tests.

--- Updates from SnPM13.0.12 ---
* snpm_pi_OneSampT
Fix: Added new random mode for sampling of sign flips (Above 53 scans the binary-coding fails as Matlab's double significand runs out of precision).

--- Updates from SnPM13.0.11 ---
* snpm_ui
Fix: When relative thresholding is set, then global calculation (user defined or mean) must be computed.

* snpm_pi_TwoSampPairT, snpm_pi_TwoSampT
Fix: Unbalanced two-sample tests were wrongly aborted.

* test_multisubpaired2cond, test_oneSample, test_twoSample, test_twosample_twocond
Updates affecting tests only.

--- Updates from SnPM13.0.10 ---
* snpm_ui
Fix: ANCOVA normalisation (with an image presenting negative mean using spm_globals).

--- Updates from SnPM13.09 ---
* snpm_cp
Increase tolerance in checking that 'Locs_vox' is integer.

--- Updates from SnPM13.08 ---
* snpm_cp, spm_append_96
Display warning message if SnPM_ST file becomes very large. Display more meaningful error message when Locs_vox is not interger.

* snpm_pp
Display warning message if SnPM_ST can not be loaded.

--- Updates from SnPM13.07 ---
* snpm
For backward-compatibility, open SPM batch window and display where to find SnPM tools when 'snpm' in called at the Matlab prompt.

* snpm_pp
Fix: slow cluster inference with variance smoothing.

* generic_test_snpm, test_oneSample, test_twoSample
Updates affecting tests only.

--- Updates from SnPM13.06 ---
* snpm_ui
Apply user-defined scaling.

* generic_test_snpm, test_oneSample, test_onesub_twocondrepl, test_twoSample
Updates affecting tests only.

--- Updates from SnPM13.05 ---
* snpm_STcalc
Cluster labels from spm_max are not necessarily continuous.

* snpm_pi_ANOVAbtwnS, snpm_pi_ANOVAwithinS, snpm_pi_Corr, snpm_pi_Corr1S, snpm_pi_OneSampT, snpm_pi_PairT, snpm_pi_TwoSampPairT, snpm_pi_TwoSampT, snpm_pi_TwoSampTss
Fix: if number of requested permutations is equal to the maximum number of possible permutations.

* snpm_pi_TwoSampTss, snpm_pi_TwoSampPairT, snpm_pi_TwoSampT
Compute possible permutations (PiCond) in a more efficient way to avoid potential memory faults for large groups.

--- Updates from SnPM13.04 ---
* snpm_ui
Apply absolute and relative masking.

--- Updates from SnPM13.03 ---
* snpm_pi_TwoSampT
Compute possible permutations (PiCond) in a more efficient way to avoid potential memory faults for large groups.

--- Updates from SnPM13.02 ---
* snpm_pp
Fix: if the extension is missing in the provided result file name, use .nii by default.
Use NaNs in thresholded map background instead of zeros.

* [deleted spm_SnPM13]
SnPM is now only accessible through the Batch window (SPM -> Tools -> SnPM).

* snpm_pi_ANOVAbtwnS, snpm_pi_ANOVAwithinS, snpm_pi_Corr1S, snpm_pi_OneSampT, snpm_pi_PairT, snpm_pi_TwoSampPairT, snpm_pi_TwoSampT
Minor update in the display.

* test/
Miror corrections affecting tests procedure only.

--- Updates from SnPM13.01 ---
* snpm_cp 
Synopsis: When half the permutations are computed (i.e. bhPerms is true) and StartPerm is 2 then the initial number of permutation could be 2 depending on sign of T0.

* snpm_test_config, snpm_test_gound_truth and generic_test_snpm
Miror corrections affecting tests procedure only.
