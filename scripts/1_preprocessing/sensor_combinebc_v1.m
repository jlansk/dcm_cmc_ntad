%% Combine gradiometers for sensor analysis
% This script combines planar gradiometers for the sensor-level analyses

subjects=BLsubs;%subjects=TWsubs
ana_dir= E.raw; %'/imaging/rowe/users/jl01/meg/dcm_cmc/meg_data/';
task='mmn';
session = 'BL';
lfile= 'ffmraeMaffffdtsss.mat';
combplanar=1;
rebasel=1;

%% Take file from NTAD preprocessing, combine planar and rebaseline
if combplanar
   for sub=1:length(subjects)
        if  ~exist((strcat(ana_dir, subjects{sub}, '/', session, '/', task, '/P', lfile)), 'file')
             ses_dir = fullfile(ana_dir,subjects{sub},session,task); % 
             S = [];
             S.D = spm_eeg_load([ses_dir '/' lfile]);
             S.mode = 'replace';
             S.prefix = 'CP';

             D = spm_eeg_combineplanar(S)
        end
   end
    lfile=['CP' lfile];
end

if rebasel
   for ss=1:length(subjects)
       S = []; S.D = spm_eeg_load([ana_dir subjects{ss} '/' session '/' task '/' lfile]);
       S.timewin = [-100 0];
       S.save = 1;
       S.prefix='b';
       D = spm_eeg_bc(S);
   end
   lfile=['b' lfile];
end