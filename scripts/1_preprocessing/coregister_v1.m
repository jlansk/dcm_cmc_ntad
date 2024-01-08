%% Automated Source Reconstruction
% JL 2023, inherited from Ece K, inherited from Alex C who inherited from Rik H

% This script coregisters each subjects T1-weighted MRI to their MEG files
% for DCM source inversion

clear all

%% switches
E = cmc_environment; %JL
filebegin='';
fileend='ffmraeMaffffdtsss.mat';
ses='AF';
bwd=E.raw;
ana = ['/' ses '/mmn/'];

if ses == 'BL' % JL
    load([E.scr filesep 'BLsubs']);
    subs=BLsubs;
elseif ses == 'AF'
    load([E.scr filesep 'AFsubs'])
    subs=AFsubs;
end

%% Parameters
epo = [-100 400]; % inversion window in ms, which is your epoch
f_win=[0.1 45];% frequency range to invert 0:all [0 4]:delta [4 8]:theta [8 12]:alpha [12 30]:beta  [30 100]:gamma
inv_typ = 'COH' % inversion method 'IID' 'GS' 'MSP'
val=1; % COH in 2, EBB in 1, IID in 3, ARD in 4, MSP in 5 , LOR in 6, GS in 7- arbitrary inversion identifier - JL change for second inversion type...
hanning=0; % hanning window switch
mesh_size = 2;   %NB: Amir has 3 ... [1-3] for coarse-medium-fine
trialtype='evoked'; %depends if it's single trial or average
inv_mods={'MEGPLANAR'};%JL - READ: make sure ece's code is commented (unless you want all modalities available to be used). Inverted fname will be first cell only modality to be inverted (NB they'll be fused i think). Leave blank to use all modalities available
inv_trls = {'DEV','REP1','REP2','REP3','REP4','REP5','REP6','REP7','REP8','REP9','REP10'};
use_headshape = 1; % switch for inclusion of digitisation points
con_win={[100 300]};%{}; % focused window,[0 500]

% Processing options
dogetmesh = 1; % process mr
docoregister = 1; % coregistration
autocoregister = 0; % use coregister if you do not have fiducial mat file
doformod = 1; % forward modelling
doinversion= 1; % invert

for ss = 1:length(subs)
    
    sub_dir=[bwd subs{ss} ana];
    D = spm_eeg_load([sub_dir filebegin fileend]);
    
    % =====================================================================
    %% Get cortical mesh
    % =====================================================================
    
    submr=[bwd subs{ss} '/sMRI/'];
    if exist(submr, 'dir') && ~isempty(submr) %EK
        l= dir(fullfile(submr, 'sMR*.nii')); %JL
        D.val = val;
        D.inv{val}.date    = strvcat(date,datestr(now,15));
        D.inv{val}.comment = {sprintf('%s_%s',inv_typ,char(inv_mods)')};    % remember inversion type and sensor configuration
        D.inv{val}.mesh    = [];
        D.inv{val}.datareg = [];
        D.inv{val}.forward = [];
        D.inv{val}.mesh.sMRI = [submr l.name];
    end
    
    if isempty(D.inv{val}.mesh.sMRI)
        error('MRI not found')
    end
    
    if exist(submr, 'dir') && ~isempty(submr)
        D.inv{val}.mesh.def  = spm_select('FPList',[submr l.name],'^y_sM.*\.nii$');% Check whether inverse-normalised surfaces present
    else
        D.inv{val}.mesh.def  = spm_select('FPList','/imaging/projects/cbu/ntad/scripts/canonical_brain/single_subj_T1.nii','^y_sM.*\.nii$');
    end
    
    if isempty(D.inv{val}.mesh.def)
        D.inv{val}.mesh = spm_eeg_inv_mesh(D.inv{val}.mesh.sMRI, mesh_size);
    end
    
    % =====================================================================
    %% Co-registration
    % =====================================================================
    newmrifid           = [];
    if autocoregister
        newmrifid.fid.pnt   = D.inv{val}.mesh.fid.fid.pnt(1:3,:); % EK force selection
    else
        fiducials_num_array=load([bwd subs{ss} '/sMRI/fiducials_num_array.mat']);
        newmrifid.fid.pnt = fiducials_num_array.fiducials_num_array;
    end
    
    newmrifid.fid.label = {'Nasion';'LPA';'RPA'};
    newmrifid.pnt       = D.inv{val}.mesh.fid.pnt;
    
    meegfid = D.fiducials;
    meegfid.pnt(find(meegfid.pnt(:,2)>0 & meegfid.pnt(:,3)<0),:) = [];
    
    D = spm_eeg_inv_datareg_ui(D, val, meegfid, newmrifid,use_headshape);
    
    
    % =====================================================================
    %% Forward Modelling
    % =====================================================================
    
    
    fprintf(1, 'Creating forward model\n');
    D.inv{val}.forward = struct([]);
    
    for ind = 1:length(D.inv{val}.datareg) %EK
        if strcmp(D.inv{val}.datareg(ind).modality,'MEG')
            D.inv{val}.forward(ind).voltype = 'Single Shell';
        elseif strcmp(D.inv{val}.datareg(ind).modality,'EEG')
            D.inv{val}.forward(ind).voltype = 'EEG BEM';
        end
    end
    
    fprintf(1, 'Computing leadfield\n');
    D = spm_eeg_inv_forward(D,val);
    
    
    % =====================================================================
    %% Inversion
    % =====================================================================
   
    D.inv{val}.inverse = [];
    D.inv{val}.inverse.trials = inv_trls ;
    D.inv{val}.inverse.type   = inv_typ;
    D.inv{val}.inverse.lpf    = f_win(1);
    D.inv{val}.inverse.hpf    = f_win(2);
    D.inv{val}.inverse.woi    = epo;% as much as the memory can handle
    D.inv{val}.inverse.Han    = hanning; %no hanning
    D.inv{val}.inverse.modality = inv_mods;
    D.inv{val}.inverse.dplot = 0;
    
    %D = spm_eeg_invert(D,val);
    
    D.save
    
end
