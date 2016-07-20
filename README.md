### SnPM: Statistical nonParametric Mapping

The <b>S</b>tatistical <b>n</b>on-<b>P</b>arametric <b>M</b>apping (SnPM) toolbox provides an extensible framework for voxel level non-parametric permutation/randomisation tests of functional Neuroimaging experiments with independent observations. 

The SnPM toolbox provides an alternative to the Statistics section of [SPM](http://www.fil.ion.ucl.ac.uk/spm/). SnPM uses the General Linear Model to construct pseudo t-statistic images, which are then assessed for significance using a standard non-parametric multiple comparisons procedure based on randomisation/permutation testing. It is most suitable for single subject PET/SPECT analyses, or designs with low degrees of freedom available for variance estimation. In these situations the freedom to use weighted locally pooled variance estimates, or variance smoothing, makes the non-parametric approach considerably more powerful than conventional parametric approaches, as are implemented in SPM. Further, the non-parametric approach is always valid, given only minimal assumptions.

#### Testing

##### Initial set up 
The first time you run the tests, you will first have to create a set of ground truth data and a file named `snpm_test_config.m` containing the reference SnPM version you used to compute the ground truth and the path to the ground truth folder. For example:
```
global testDataDir;
testDataDir = '~/snpm_test_data';
global SnPMrefVersion;
SnPMrefVersion = 'SnPM8';
```

##### Run the tests
Then, the tests can be started (from the test data folder) with:
```
import matlab.unittest.TestSuite;
suite = TestSuite.fromFolder(fullfile(spm_str_manip(which('snpm'), 'h'), 'test'));
result = run(suite);
```