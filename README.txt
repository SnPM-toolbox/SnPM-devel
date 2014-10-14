SnPM13: Bugs & Fixes

This file describes the bugs that have been reported, along with the appropriate fixes. 
Updated versions of appropriate SnPM functions are available from in the snpm13_updates: 
http://warwick.ac.uk/snpm/distribution/snpm13_updates

--- SnPM13.1.01 ---
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