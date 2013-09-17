% Compare results (t-map and positive effects filtered map T>3) of
% interactive and batch one subject, 2 conditions with replication
function test_suite = test_onesub_twocondrepl()
    global TEST;
    TEST = true;

    global parentDataDir;
    parentDataDir = 'M:\Data\snpm_test_data';
    global testDataDir;
    testDataDir = fullfile(parentDataDir, 'data', 'PET_motor');

    initTestSuite;
    
    clear global TEST;
end

% No covariate, no variance smoothing
function test_onesub_twocondrepl_1()
    generic_for_all();
end
function matlabbatch = onesub_twocondrepl_1_batch(matlabbatch)
end

% With variance smoothing
function test_onesub_twocondrepl_var()
    generic_for_all();
end
function matlabbatch = onesub_twocondrepl_var_batch(matlabbatch)
    matlabbatch{1}.cfg_snpm.Design.TwoSampTss.vFWHM = [12 12 12];
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
    matlabbatch{1}.cfg_snpm.Design.TwoSampTss.dir = {batchResDir};
    matlabbatch{1}.cfg_snpm.Design.TwoSampTss.P = {
         fullfile(testDataDir, 's8np01160em01R.img,1')
         fullfile(testDataDir, 's8np01160em02R.img,1')
         fullfile(testDataDir, 's8np01160em03R.img,1')
         fullfile(testDataDir, 's8np01160em04R.img,1')
         fullfile(testDataDir, 's8np01160em05R.img,1')
         fullfile(testDataDir, 's8np01160em06R.img,1')
         fullfile(testDataDir, 's8np01160em07R.img,1')
         fullfile(testDataDir, 's8np01160em08R.img,1')
         fullfile(testDataDir, 's8np01160em09R.img,1')
         fullfile(testDataDir, 's8np01160em10R.img,1')
         fullfile(testDataDir, 's8np01160em11R.img,1')
         fullfile(testDataDir, 's8np01160em12R.img,1')
         };
    matlabbatch{1}.cfg_snpm.Design.TwoSampTss.condidx = [1 2 1 2 1 2 1 2 1 2 1 2];
    matlabbatch{1}.cfg_snpm.Design.TwoSampTss.Tss_repc = 6;
    matlabbatch{1}.cfg_snpm.Design.TwoSampTss.TwosampTss_Block = 4;

    matlabbatch{1}.cfg_snpm.Design.OneSampT.dir = {batchResDir};
    matlabbatch{1}.cfg_snpm.Design.OneSampT.P = {
                                                 fullfile(testDataDir, 'con_0001.img,1')
                                                 fullfile(testDataDir, 'con_0001.img,1')
                                                 fullfile(testDataDir, 'su_control03', 'cn_sess1', 'con_0001.img,1')
                                                 fullfile(testDataDir, 'su_control04', 'cn_sess1', 'con_0001.img,1')
                                                 fullfile(testDataDir, 'su_control05', 'cn_sess1', 'con_0001.img,1')
                                                 };
    
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