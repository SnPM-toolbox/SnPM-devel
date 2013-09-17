% Compare results (t-map and positive effects filtered map T>3) of
% interactive and batch one-sample t-tests

function test_suite = test_oneSample()
    global TEST;
    TEST = true;

    global parentDataDir;
    parentDataDir = 'M:\Data\snpm_test_data';
    global testDataDir;
    testDataDir = fullfile(parentDataDir, 'data');

    initTestSuite;
    
    clear global TEST;
end

% No covariate, no variance smoothing
function test_onesample_1()
    generic_for_all();
end
function matlabbatch = onesample_1_batch(matlabbatch)
end

% With 1 covariate
function test_onesample_cov()
    generic_for_all();
end
function matlabbatch = onesample_cov_batch(matlabbatch)
    %matlabbatch{1}.cfg_snpm.Design.OneSampT.covariate.cov_Val = [1 5 2 21 0];
    matlabbatch{1}.cfg_snpm.Design.OneSampT.mcov.c = [1 5 2 21 0];
    matlabbatch{1}.cfg_snpm.Design.OneSampT.mcov.cname = 'age';
end

% With 3 covariates
function test_onesample_cov3()
    generic_for_all();
end
function matlabbatch = onesample_cov3_batch(matlabbatch)
    matlabbatch{1}.cfg_snpm.Design.OneSampT.mcov(1).c = [1 1 2 3 1];
    matlabbatch{1}.cfg_snpm.Design.OneSampT.mcov(1).cname = 'age';
    matlabbatch{1}.cfg_snpm.Design.OneSampT.mcov(2).c = [0 21 15 18 3];
    matlabbatch{1}.cfg_snpm.Design.OneSampT.mcov(2).cname = 'height';
    matlabbatch{1}.cfg_snpm.Design.OneSampT.mcov(3).c = [-1 -0.5 -1 1 0];
    matlabbatch{1}.cfg_snpm.Design.OneSampT.mcov(3).cname = 'width';
end

% With variance smoothing
function test_onesample_var()
    generic_for_all();
end
function matlabbatch = onesample_var_batch(matlabbatch)
    matlabbatch{1}.cfg_snpm.Design.OneSampT.vFWHM = [6 6 6];
end

% With approximate test
function test_onesample_approx()
    rand('seed',200);
    generic_for_all();
end
function matlabbatch = onesample_approx_batch(matlabbatch)
    global testDataDir;
    matlabbatch{1}.cfg_snpm.Design.OneSampT.P(end+1:end+8) = {
         fullfile(testDataDir, 'su_control06', 'cn_sess1', 'con_0001.img,1')
         fullfile(testDataDir, 'su_control07', 'cn_sess1', 'con_0001.img,1')
         fullfile(testDataDir, 'su_control08', 'cn_sess1', 'con_0001.img,1')
         fullfile(testDataDir, 'su_control09', 'cn_sess1', 'con_0001.img,1')
         fullfile(testDataDir, 'su_control10', 'cn_sess1', 'con_0001.img,1')
         fullfile(testDataDir, 'su_control11', 'cn_sess1', 'con_0001.img,1')
         fullfile(testDataDir, 'su_control12', 'cn_sess1', 'con_0001.img,1')
         fullfile(testDataDir, 'su_control13', 'cn_sess1', 'con_0001.img,1')
         };
    matlabbatch{1}.cfg_snpm.Design.OneSampT.nPerm = 100;
end

function assert_check(batchResDir, interResDir, interNewResDir)
    % Compare t-maps and filtered maps
    batch_tmap = spm_select('FPList', batchResDir, '^snpmT\+\.img');
    batch_filtmap = spm_select('FPList', batchResDir, '^SnPM_filtered_10none.*\.nii');
    
    inter_tmap = spm_select('FPList', interResDir, '^snpmT\+\.img');
    inter_filtmap = spm_select('FPList', interResDir, '^SnPMt_filtered_10none\.img');
    
    assertEqual(spm_read_vols(spm_vol(batch_tmap)), spm_read_vols(spm_vol(inter_tmap)))
    assertEqual(spm_read_vols(spm_vol(batch_filtmap)), spm_read_vols(spm_vol(inter_filtmap)))
    
    % Useless, do not keep spm2-like interface
%     % Equality of new interactive with "ground truth" from previous snpm
%     % version
%     inter_new_tmap = spm_select('FPList', interNewResDir, '^snpmT\+\.img');
%     inter_new_filtmap = spm_select('FPList', interNewResDir, '^SnPMt_filtered_10none\.img');
%     
%     assertEqual(spm_read_vols(spm_vol(inter_new_tmap)), spm_read_vols(spm_vol(inter_tmap)))
%     assertEqual(spm_read_vols(spm_vol(inter_new_filtmap)), spm_read_vols(spm_vol(inter_filtmap)))
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
    matlabbatch{1}.cfg_snpm.Design.OneSampT.dir = {batchResDir};
    matlabbatch{1}.cfg_snpm.Design.OneSampT.P = {
                                                 fullfile(testDataDir, 'su_control01', 'cn_sess1', 'con_0001.img,1')
                                                 fullfile(testDataDir, 'su_control02', 'cn_sess1', 'con_0001.img,1')
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