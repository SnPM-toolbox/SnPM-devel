import matlab.unittest.TestSuite;
suite = TestSuite.fromFolder(fullfile(spm_str_manip(which('snpm'), 'h'), 'test'));
result = run(suite);