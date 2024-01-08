function DCM = dcm_cmc_gen_v1(DCMa, raw, ana, sub, rawlast)

% Generating the alternate model structures
%==========================================================================

% DCMa is your model specification of A B and C matrices. 
% raw is the start of your file address. 
% sub is subject name
% rawlast is end of file address of data file to be inverted
% ana is your analysis directory

        
%--------------------------------------------------------------------------
%% ---- Set A and C matrices of DCM from DCMa       
DCM = DCMa; %specifies A B and C matrices

% Data filename
%--------------------------------------------------------------------------

DCM.xY.Dfile =  [raw filesep sub rawlast];


%% ----DCM model parameters ----
% Parameters and options used for setting up model
%--------------------------------------------------------------------------
DCM.options.analysis = 'ERP'; % analyze evoked responses
DCM.options.model    = 'CMC' %'CMM_NMDA'; % CMC, TFM
DCM.options.Tdcm(1)  = 0;     % start of peri-stimulus time to be modelled
DCM.options.Tdcm(2)  = 300;   % end of peri-stimulus time to be modelled
DCM.options.Nmodes   = 8;     % number of modes for data selection
DCM.options.h        = 1;     % number of DCT components
DCM.options.onset    = 60;    % onset
DCM.options.dur     = 16;     % duration
DCM.options.D        = 1;     % downsampling
DCM.options.trials   = [1 2 3 4 5 6]; % index of trials you want to use

DCM.xY.modality='MEGPLANAR';
DCM.options.spatial  = 'ECD'; %LFP 

DCM.options.location = 1;     % optimising source location
DCM.options.han      = 1;     % applying hanning window
DCM.options.symmetry=1;


% Location priors for dipoles
%--------------------------------------------------------------------------
% ORIGINAL from T Cope et al 2022: DCM.Lpos = [[-42; -22; 7] [ -61; -32; 8] [ -46; 20; 8] [ -49; -38; 38] [ 46; -14; 8] [ 59; -25; 8] [ 46; 20; 8] [ 57; -38; 42]];

%coords from Lappe et al. 2013 for IPC coordinates
DCM.Lpos = [[-42; -22; 7] [ -61; -32; 8] [ -46; 20; 8] [-58; -27; 30] [ 46; -14; 8] [ 59; -25; 8] [ 46; 20; 8] [59; -41; 30]]; 
DCM.Sname  = {'left A1'; 'left STG'; 'left IFG'; 'left IPC'; 'right A1'; 'right STG'; 'right IFG'; 'right IPC'};
Nareas    = size(DCM.Lpos,2);



% Specify B matrix
%--------------------------------------------------------------------------

DCM.B{1} = DCM.A{1} + DCM.A{2} + DCM.A{3};  % model specification
DCM.B{1}=logical(DCM.B{1});
DCM.B{2}=logical(DCM.B{1});

        
% Between trial effects
%--------------------------------------------------------------------------
% Trials defined as deviant standard
%% ---- Condition-specific (trial) effects ---- from Garrido et al. 2009

DCM.xU.X(:,1) = [1; 0.367879; 0.135335; 0.049787; 0.018316; 0.006738]
DCM.xU.name{1} = ('1-6_expo');
DCM.xU.X(:,2)= [0; 0.541341; 0.146525; 0.029750; 0.005374; 0.000908] 
DCM.xU.name{2}={'1-6_phas'};

DCM.xU.dt = 0.002; 
        
% Define priors
%--------------------------------------------------------------------------

[pE,pC]  = spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,DCM.options.model);
DCM.M.pE = pE;
DCM.M.pC = pC;

DCM.name = [ana filesep 'DCM_' sub '_' DCMa.name];
       

end
