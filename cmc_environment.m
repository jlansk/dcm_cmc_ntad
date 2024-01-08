function E = cmc_environment()

%% Set up environment
E.Fspm = '/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7771/';
addpath(E.Fspm);
spm('defaults', 'EEG');

E.scr='/imaging/rowe/users/jl01/meg/dcm_cmc/scripts/to_publish/';
addpath(genpath(E.scr));

E.raw='/imaging/rowe/users/jl01/meg/dcm_cmc/meg_data/';

E.anaB =   '/imaging/rowe/users/jl01/meg/dcm_cmc/Oct21_0300_BL_lappeIPC_he6/64steps/full/';
E.anaL =  '/imaging/rowe/users/jl01/meg/dcm_cmc/Oct21_0300_AF_lappeIPC_he6/64steps/full/'

end