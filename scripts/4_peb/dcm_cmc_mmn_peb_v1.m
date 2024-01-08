%% run between group PEB of winning model - Baseline comparison
%  1. Compare baseline group differences (controls vs patients)
%  2. Compare longitudinal group differences (patients; baseline vs
%  follow-up)

%% Switches
clear all
pebBg = 1;
pebLg = 1; dayscov = 1;
agecov=1;
savepebfig=1;

%% Set up environment and variables
E = cmc_environment;
scr=E.scr; 
anaB = E.anaB;
anaL=E.anaL;
anaF =  [scr filesep 'figures'];

Bbmr = load([anaB filesep 'RCM_BMC_BMA_may22']); %or BMC_BMA_2 and rcm.mat
maxind = 6;
BGCM=spm_dcm_load(Bbmr.RCM(:,maxind));
Lbmr=load([anaL filesep 'RCM_BMC_BMA_BLAFonlysubs.mat']);
LGCM=spm_dcm_load(Lbmr.RCM(:,maxind));

%% Baseline Group comparison
% specify PEB model
for f=1:length(BGCM)
    Bfiles{f}=BGCM{f}.name;
end
Bsubs=extractBetween(Bfiles, 'DCM_', '_');

if pebBg == 1
    
    X=zeros(length(Bsubs),2);
    X(:,1)=ones; %first collumn to ones to model overall mean
    
    con = find(contains(Bsubs, 'C'));
    pat = find(contains(Bsubs, 'P'));
    
    X(1:con(end),2)=0; %set controls (index 1:14) to 0 ;
    X_labels={'mean','cons->pats'};
    X(pat:end,2)=1;
    
    M=struct(); %a structure to specify the PEB settings
    M.Q = 'all'; % between-subject variability will be estimated for each connection
    M.X=X; %design matrix
    M.Xnames = X_labels;
    
    field{1}={'G(:,1)'}; %superficial pyramidal cell gain
    field{2}={'A{2}','A{3}'}; % extrinsic connectivity between pyramidal cells
    
    for f = 1:length(field)
        [PEB{f}, RCM1{f}] = spm_dcm_peb(BGCM, M, field{f}); % run the PEB
        
        %post-hoc DCM of hemispherically-symetrical model space
       
        np = length(PEB{f}.Pnames);
        nr = size(BGCM{1}.Lpos,2)
        Mall{f}=spm_perm_mtx(np); Mall{f}=full(Mall{f});
        
        count=0;
        if np==nr
            for ss=1:size(Mall{f},1)
                LHS{f} = [Mall{f}(ss,1:4)];
                RHS{f} = [Mall{f}(ss,5:8)];
                if all(LHS{f} == RHS{f})
                    count=count+1;
                    sym_idx{f}(count) = ss;
                end
            end
        elseif np==nr*2
            for ss=1:size(Mall{f},1)
                LHS{f} = [Mall{f}(ss,1:4), Mall{f}(ss,9:12)];
                RHS{f} = [Mall{f}(ss,5:8), Mall{f}(ss,13:16)];
                if all(LHS{f} == RHS{f})
                    count=count+1;
                    sym_idx{f}(count) = ss;
                end
            end
        end
        
        Msymm{f}=Mall{f}(sym_idx{f},:);
        RMAmodel{f} = spm_dcm_peb_bmc(PEB{f},  Msymm{f});   %
        
    end
    
    %% save figs
    % spm_dcm_peb_review(RMAmodel{f}, BGCM');
    if savepebfig
        
        for f =1:2   
            spm_dcm_peb_review_fig_jl(RMAmodel{f}, 2, 0.95, 2)
             
            if f==1
                xticklabels({'left A1', 'left STG', 'left IFG', 'left IPC', 'right A1', 'right STG', 'right IFG', 'right IPC'});
            elseif f==2
                for x=1:length(RMAmodel{f}.Pnames)
                    xticknames{1, x}=pname_to_string(RMAmodel{f}.Pnames{x}, BGCM{1}.Sname);
                    xticknames{2, x}=strcat(extractBefore(xticknames{1,x}, '-'), ':  ', extractBetween(xticknames{1,x}, 'matrix  ', 'to'), 'to', xticknames{1,x}(end-3:end));
                    xticknames{3, x}=strcat(extractBetween(xticknames{1,x}, 'matrix  ', 'to'), 'to', xticknames{1,x}(end-3:end));
                    
                end
                xticklabels([xticknames{3,:}]);
            end
            
            exportgraphics(gcf, [anaF '/group_' labels{f} '_95.png'], 'Resolution', '720');
        end
    end
end

%% Longitudinal PEB analysis

if pebLg == 1
  
    load([anaL '/Lsubs.mat']);
    LX=zeros(2*length(Lsubs),2);
    LX(:,1)=ones; %first collumn to ones to model overall mean
    
    LX(1:length(Lsubs),2)=zeros;
    
    if dayscov == 1
        load([anaL '/days_af-bl.mat']) %should load into a variable called 'days'
        
        
        for ss=1:length(Lsubs)
            idx=find(contains(days(:,1),Lsubs(ss)));
            LX(length(Lsubs)+ss,3) = days{idx,2}/365.2425;
        end
        LX(length(Lsubs)+1:length(Lsubs)*2,3) = zscore(LX(length(Lsubs)+1:length(Lsubs)*2,3));
        LX(length(Lsubs)+1:end,2) = 1;
        LX_labels={'mean','BL->AF', 'days'};

    else
        LX(length(Lsubs)+1:end,2) = 1;
        LX_labels={'mean','BL->AF'};

    end
    
    
    %PEB settings - PZ:"there are 2 prerequisites for a PEB analysis: a
    %between-subject design matrix and the collated GCM file"
    
    LM=struct(); %a structure to specify the PEB settings
    LM.Q = 'all'; % between-subject variability will be estimated for each connection
    LM.X=LX; %design matrix
    LM.Xnames = LX_labels;
    
    [LPEB{1}, RCMG] = spm_dcm_peb(LGCM, LM, {'G(:,1)'});
    [LPEB{2}, RCMA] = spm_dcm_peb(LGCM, LM, {'A{2}', 'A{3}'});
    
    %% Model longitudinal connections >.95 in baseline model comparison
    MG = [1,0,1,0,1,0,1,0;1,0,0,0,1,0,0,0;0,0,1,0,0,0,1,0;0,0,0,0,0,0,0,0]; %symmetric winning parameters from BL 
    LRMAmodelG = spm_dcm_peb_bmc(LPEB{1},  MG);   %
    %spm_dcm_peb_review(LRMAmodelG, LGCM');
    
    MA = [0,1,1,1,0,1,1,1,0,0,1,0,0,0,1,0;0,1,1,0,0,1,1,0,0,0,1,0,0,0,1,0;0,1,0,1,0,1,0,1,0,0,1,0,0,0,1,0;0,1,0,0,0,1,0,0,0,0,1,0,0,0,1,0;0,0,1,1,0,0,1,1,0,0,1,0,0,0,1,0;0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0;0,0,0,1,0,0,0,1,0,0,1,0,0,0,1,0;0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0;0,1,1,1,0,1,1,1,0,0,0,0,0,0,0,0;0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0;0,1,0,1,0,1,0,1,0,0,0,0,0,0,0,0;0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0;0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,0;0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0;0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
    LRMAmodelA = spm_dcm_peb_bmc(LPEB{2},  MA);   %
    %spm_dcm_peb_review(LRMAmodelA, LGCM');
    
    spm_dcm_peb_review_fig_jl(LRMAmodelG, 2, 0.95, 2); xticklabels({'left A1', 'left STG', 'left IFG', 'left IPC', 'right A1', 'right STG', 'right IFG', 'right IPC'});
   
    spm_dcm_peb_review_fig_jl(LRMAmodelA, 2, 0.95, 2);
    for x=1:length(LRMAmodelA.Pnames)
        xticknames{1, x}=pname_to_string(LRMAmodelA.Pnames{x}, BGCM{1}.Sname);
        xticknames{2, x}=strcat(extractBetween(xticknames{1,x}, 'matrix  ', 'to'), 'to', xticknames{1,x}(end-3:end));
        %xticknames{2, x}=strcat(extractBefore(xticknames{1,x}, '-'), ':  ', extractBetween(xticknames{1,x}, 'matrix  ', 'to'), 'to', xticknames{1,x}(end-3:end));
    end
    xticklabels([xticknames{2,:}])
    
end

%% repetition - supplementary analysis

%% Modulatory intrinsic

PEB_N = spm_dcm_peb(BGCM, M, {'N{2}'});
RMA_N = spm_dcm_peb_bmc(PEB_N, Msymm{1})

LPEB_N = spm_dcm_peb(LGCM, LM, {'N{2}'});
LRMA_N = spm_dcm_peb_bmc(LPEB_N, [0,0,1,0,0,0,1,0;0,0,0,1,0,0,0,1;0,0,1,1,0,0,1,1;0,0,0,0,0,0,0,0]);

spm_dcm_peb_review_fig_jl(RMA_N, 2, 0.95, 2);
xticklabels({'left A1', 'left STG', 'left IFG', 'left IPC', 'right A1', 'right STG', 'right IFG', 'right IPC'})
exportgraphics(gcf, [anaF '/group_N_95_symm.png'], 'Resolution', '720');

spm_dcm_peb_review_fig_jl(LRMA_N, 2, 0.95, 2);
xticklabels({'left A1', 'left STG', 'left IFG', 'left IPC', 'right A1', 'right STG', 'right IFG', 'right IPC'})
exportgraphics(gcf, [anaF '/L_N_95_symm.png'], 'Resolution', '720');

%% Modulatory extrinsic
%PEB_B = spm_dcm_peb(BGCM, M, {'B{1}'});%
PEB_B = spm_dcm_peb(BGCM, M, {'B{1}(2,1)','B{1}(3,2)','B{1}(4,2)','B{1}(6,5)','B{1}(7,6)','B{1}(8,6)','B{1}(3,4)','B{1}(7,8)'});%'B{1}(1,2)','B{1}(2,3)','B{1}(2,4)','B{1}(5,6)','B{1}(6,7)','B{1}(6,8)','B{1}(4, 3)','B{1}(8,7)'});
RMA_B = spm_dcm_peb_bmc(PEB_B, Msymm{1})

spm_dcm_peb_review_fig_jl(RMA_B, 2, 0.95, 2);
xticknames=xticklabels_jl_peb(PEB_B, BGCM{1});
xticklabels([xticknames{2,:}])
exportgraphics(gcf, [anaF '/group_B1f_95_symm.png'], 'Resolution', '720');

% LPEB_B = spm_dcm_peb(LGCM, LM, {'B{1}'});
LPEB_B = spm_dcm_peb(LGCM, LM, {'B{1}(2,1)','B{1}(3,2)','B{1}(4,2)','B{1}(6,5)','B{1}(7,6)','B{1}(8,6)','B{1}(3,4)','B{1}(7,8)'});%'B{1}(1,2)','B{1}(2,3)','B{1}(2,4)','B{1}(5,6)','B{1}(6,7)','B{1}(6,8)','B{1}(4, 3)','B{1}(8,7)'});
LRMA_B = spm_dcm_peb_bmc(LPEB_B, [0,0, 1, 1, 0, 0, 1, 1; 0, 0, 1, 0, 0, 0, 1, 0; 0,0,0,1,0,0,0,1; 0,0,0,0,0,0,0,0]);%,0,0,0,1,0,0;0,0,0,0,0,0,0,0])

spm_dcm_peb_review_fig_jl(LRMA_B, 2, 0.95, 2);
xticklabels([xticknames{2,:}])
exportgraphics(gcf, [anaF '/L_B1f_95_symm.png'], 'Resolution', '720');

%% is A or B & G or N most important for explaining effects

F_AB(1)=PEB{2}.F; F_AB(2)=PEB_B.F;
plot_bmc_F_jl(F_AB, {'A', 'B'}, ('Group extrinsic BMC'));
exportgraphics(gcf, [anaF '/group_AB_BMC.png'], 'Resolution', '720');


F_GN(1)=PEB{1}.F; F_GN(2)=PEB_N.F
plot_bmc_F_jl(F_GN, {'G', 'N'}, ('Group intrinsic BMC'));
exportgraphics(gcf, [anaF '/group_GN_BMC.png'], 'Resolution', '720');


LF_AB(1)=LPEB{2}.F; LF_AB(2)=LPEB_B.F;
plot_bmc_F_jl(LF_AB, {'A', 'B'}, ('Longitudinal extrinsic BMC'));
exportgraphics(gcf, [anaF '/L_AB_BMC.png'], 'Resolution', '720');

LF_GN(1)=LPEB{1}.F; LF_GN(2)=LPEB_N.F;
plot_bmc_F_jl(LF_GN, {'G', 'N'}, ('Longitudinal intrinsic BMC'));
exportgraphics(gcf, [anaF '/L_GN_BMC.png'], 'Resolution', '720');


%% Age as covariate for group analysis
if agecov==1
    
    load([anaB '/BLageMEG.mat']);
    for ss=1:length(Bsubs)
        idx=find(contains(BLage(:,1),Bsubs(ss)));
        Xage(ss,1) = BLage(idx,1);
        Xage(ss,2) = BLage(idx,2); % should be the value corresponding to the peb_subs
    end
    
    for ss=1:length(Bsubs)
        if strcmp(Xage{ss,1}, Bsubs{ss})
            X(ss,3)=Xage{ss,2}; %age if needed?
        else
            disp('subject orders do not match')
        end
    end
    X_labels={'mean','cons>pats','age'};
    X(:,3:end)=zscore(X(:,3:end));
    
    Mage=struct(); %a structure to specify the PEB settings
    Mage.Q = 'all'; % between-subject variability will be estimated for each connection
    Mage.X=X; %design matrix
    Mage.Xnames = X_labels;
    
    [PEB_age_G, RCM1_age_G] = spm_dcm_peb(BGCM, Mage, {'G(:,1)'});
    RMAmodel_age_G = spm_dcm_peb_bmc(PEB_age_G,  MG);%
    spm_dcm_peb_review_fig_jl(RMAmodel_age_G, 2, 0.95, 2);
    
    [PEB_age_A, RCM1_age_A] = spm_dcm_peb(BGCM, Mage, {'A{2}','A{3}'});
    RMAmodel_age_A = spm_dcm_peb_bmc(PEB_age_A,  MA);%
    spm_dcm_peb_review_fig_jl(RMAmodel_age_A, 2, 0.95, 2);
end
