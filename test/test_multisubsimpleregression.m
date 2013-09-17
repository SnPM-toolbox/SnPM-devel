function test_suite = test_multisubsimpleregression()
    global TEST;
    TEST = true;

    global parentDataDir;
    parentDataDir = 'M:\Data\snpm_test_data';
    global testDataDir;
    testDataDir = fullfile(parentDataDir, 'data');

    initTestSuite;
    
    clear global TEST;
end

% No variance smoothing, no approximate test
function test_multisubsimpleregression_1()
    generic_for_all();
end
function matlabbatch = multisubsimpleregression_1_batch(matlabbatch)
end

% With variance smoothing
function test_multisubsimpleregression_var()
    generic_for_all();
end
function matlabbatch = multisubsimpleregression_var_batch(matlabbatch)
    matlabbatch{1}.cfg_snpm.Design.Corr.vFWHM = [6 6 6];
end

% With approximate test
function test_multisubsimpleregression_approx()
    rand('seed',200);
    generic_for_all();
end
function matlabbatch = multisubsimpleregression_approx_batch(matlabbatch)
    global testDataDir;
    matlabbatch{1}.cfg_snpm.Design.Corr.P(end+1:end+8) = {
         fullfile(testDataDir, 'su_control06', 'cn_sess1', 'con_0001.img,1')
         fullfile(testDataDir, 'su_control07', 'cn_sess1', 'con_0001.img,1')
         fullfile(testDataDir, 'su_control08', 'cn_sess1', 'con_0001.img,1')
         fullfile(testDataDir, 'su_control09', 'cn_sess1', 'con_0001.img,1')
         fullfile(testDataDir, 'su_control10', 'cn_sess1', 'con_0001.img,1')
         fullfile(testDataDir, 'su_control11', 'cn_sess1', 'con_0001.img,1')
         fullfile(testDataDir, 'su_control12', 'cn_sess1', 'con_0001.img,1')
         fullfile(testDataDir, 'su_control13', 'cn_sess1', 'con_0001.img,1')
         };
    matlabbatch{1}.cfg_snpm.Design.Corr.CovInt = [1 3 5 0 2 6 7 2 1 -1 2 3 1];
    matlabbatch{1}.cfg_snpm.Design.Corr.nPerm = 100;
end

function assert_check(batchResDir, interResDir, interNewResDir)
    % Compare t-maps and filtered maps
    batch_tmap = spm_select('FPList', batchResDir, '^snpmT\+\.img');
    batch_filtmap = spm_select('FPList', batchResDir, '^SnPM_filtered_10none.*\.nii');
    
    inter_tmap = spm_select('FPList', interResDir, '^snpmT\+\.img');
    inter_filtmap = spm_select('FPList', interResDir, '^SnPMt_filtered_10none\.img');
    
    assertEqual(spm_read_vols(spm_vol(batch_tmap)), spm_read_vols(spm_vol(inter_tmap)))
    assertEqual(spm_read_vols(spm_vol(batch_filtmap)), spm_read_vols(spm_vol(inter_filtmap)))
end

% Must *not* have "test" in the name of this generic function
function generic_for_all()
    global TEST;
    TEST = true;

    global parentDataDir;
    global testDataDir;
    
    st = dbstack();
    testDirName = lower(st(2).name(6:end));
    batchResDir = fullfile(parentDataDir, 'results', 'batch', testDirName);
    matlabbatch = create_basis_matlabbatch(batchResDir,testDataDir);
    eval(['matlabbatch = ' testDirName '_batch(matlabbatch)']);
    spm_jobman('run', matlabbatch);
    
    interResDir = fullfile(spm_str_manip(batchResDir,'hh'), 'interactive', testDirName);
    interNewResDir = fullfile(spm_str_manip(batchResDir,'hh'), 'interactive_new', testDirName);
    assert_check(batchResDir, interResDir, interNewResDir);
    
    clear global TEST;
end

function matlabbatch = create_basis_matlabbatch(batchResDir,testDataDir)
    matlabbatch{1}.cfg_snpm.Design.Corr.dir = {batchResDir};
    matlabbatch{1}.cfg_snpm.Design.Corr.P = {
                                                 fullfile(testDataDir, 'su_control01', 'cn_sess1', 'con_0001.img,1')
                                                 fullfile(testDataDir, 'su_control02', 'cn_sess1', 'con_0001.img,1')
                                                 fullfile(testDataDir, 'su_control03', 'cn_sess1', 'con_0001.img,1')
                                                 fullfile(testDataDir, 'su_control04', 'cn_sess1', 'con_0001.img,1')
                                                 fullfile(testDataDir, 'su_control05', 'cn_sess1', 'con_0001.img,1')
                                                 };
    matlabbatch{1}.cfg_snpm.Design.Corr.CovInt = [1 3 5 0 2];
    
    % Compute
    matlabbatch{2}.cfg_snpm.snpm_bch_cp.snpmcfg(1) = cfg_dep;
    matlabbatch{2}.cfg_snpm.snpm_bch_cp.snpmcfg(1).tname = 'SnPMcfg.mat configuration file';
    matlabbatch{2}.cfg_snpm.snpm_bch_cp.snpmcfg(1).tgt_spec = {};
    matlabbatch{2}.cfg_snpm.snpm_bch_cp.snpmcfg(1).sname = 'MultiSub: One Sample T test on diffs/contrasts: SnPMcfg.mat configuration file';
    matlabbatch{2}.cfg_snpm.snpm_bch_cp.snpmcfg(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{2}.cfg_snpm.snpm_bch_cp.snpmcfg(1).src_output = substruct('.','SnPMcfg');
    
    % Results   
    matlabbatch{3}.cfg_snpm.Infer.SnPMmat(1) = cfg_dep;
    matlabbatch{3}.cfg_snpm.Infer.SnPMmat(1).tname = 'SnPM.mat results file';
    matlabbatch{3}.cfg_snpm.Infer.SnPMmat(1).tgt_spec = {};
    matlabbatch{3}.cfg_snpm.Infer.SnPMmat(1).sname = 'Compute: SnPM.mat results file';
    matlabbatch{3}.cfg_snpm.Infer.SnPMmat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1});
    matlabbatch{3}.cfg_snpm.Infer.SnPMmat(1).src_output = substruct('.','SnPM');
    %matlabbatch{3}.cfg_snpm.Infer.Thr.Vox.VoxSig.TFth = 3;
    matlabbatch{3}.cfg_snpm.Infer.Thr.Vox.VoxSig.Pth = 0.10;
    matlabbatch{3}.cfg_snpm.Infer.Tsign = 1;
    matlabbatch{3}.cfg_snpm.Infer.WriteFiltImg.name = 'SnPM_filtered_10none.nii';
end