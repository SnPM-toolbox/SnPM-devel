% Create dummy scans for testing purposes.

%_______________________________________________________________________
% Copyright (C) 2013-14 The University of Warwick
% Id: snpm_create_test_scans.m  SnPM8 2014/01/31
% Camille Maumet
function snpm_create_test_scans()
snpm_test_config;
global testDataDir;

if isempty(testDataDir)
    error('Test data directory not set, please update snpm_test_config');
end

buttonName = questdlg('Current test data will be overwritten. Is that ok?',...
    'Overwrite','yes','no','no');

switch buttonName,
    case 'no',
        error('Re-computation of test data aborted!');
    case 'yes',
        
        % We need 17 scans in one-sample test with approximated
        % distribution
        nScans = 17;
        gr2nScans = 3;
        
        % Create 'nScans' dummy (small) volumes to be used for testing.
        for i = (nScans+1):(nScans+gr2nScans)
            grSuffix = '';
            if i > nScans
                grSuffix = 'gr2_';
            end
            
            vol.fname    = fullfile(testDataDir, ['test_data_' grSuffix num2str(i, '%02.0f') '.nii']);
            vol.descrip  = ['Test data ' num2str(i, '%02.0f')];
            vol.mat      =  [   -2     0     0    80; ...
                0     2     0  -114; ...
                0     0     4   -56; ...
                0     0     0     1];
            vol.dim      = [10 10 10]; % Create small images
            
            vol  = spm_create_vol(vol);
            
            % Noise data...
            data = normrnd(0,1, vol.dim);
            
            if i < nScans
                % Add some big effect (for FWE detections)
                data(2:4,2:4,2:5) = normrnd(10,1, [3 3 4]);
            end
            
            vol = spm_write_vol(vol, data);
        end
end
end