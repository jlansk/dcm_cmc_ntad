%% DCM using the cmc model to model Alzheimer's disease progression
%  1. Inverts the full DCM for each individual and check fits
%  2. Runs a bayesian model reduction to compare models to identify how
%  parietal node connects

%% Set up environment and variables
clearvars
E = cmc_environment;
scr= E.scr; 
anaL= E.anaL;
anaB= E.anaB;
raw= E.raw;

load([scr '/AFsubs.mat']);
subjects=AFsubs;
rawlast = '/AF/mmn/ffmraeMaffffdtsss.mat'; % file that has had all NTAD preprocessing steps applied) and link to mri

%% switches
invert = 1;
plot1 = 1;
bmc = 1 ;
plot2 = 1;

%% Set up  model space

%create models from my model space with different A and C matrices
[DCMa] = dcm_cmc_gen_model_space_v1();

% specify rest of parameters for the full DCM model
for ss=1:length(subjects)
    dcms{ss} = dcm_cmc_gen_v1(DCMa{end}, raw, anaL, subjects{ss}, rawlast);
end

dcms=dcms';
GCM=spm_dcm_load(dcms);

%% Invert first-level DCMs
%  subject that do not converge go into the 'failed' variable

if invert == 1
    
    failed={};
    parfor ss=1:length(subjects)
        try
            spm_dcm_erp(GCM{ss});
        catch
            failed{ss}={strcat('failed_', num2str(ss))};
            continue;
        end
    end
    
    fail = find(~cellfun(@isempty,failed));
    failed_subs = subjects(fail);
    save([anaL '/failedsubs.mat'],failed_subs)
    
    for ss=1:length(fail)
        indx=fail(ss)
        spm_dcm_erp(GCM{indx})
    end
end
%% Remove failed or bad fits
% From the variable 'failed_subs' created during the inversion and from
% seeing which have NANs in the H matrix (as shown through above plots) we
% remove those failed subjects
BLbmr=load([anaB filesep 'RCM_BMC_BMA_may22']); %or BMC_BMA_2 and rcm.mat

if ~exist('maxind', 'var') %nb max ind for AF is also 6
    maxind=6; %or BMC_BMA_2 and rcm.mat
end

if ~exist('failed_subs', 'var')
    load([anaL '/failedsubs.mat']) %failed_subs={'C1012'}
end


count=0;
for ss = 1:length(subjects)
    try
        DCM =spm_dcm_load([anaL filesep 'DCM_' subjects{ss} '_full.mat']);
        DCM = DCM{1};
        if any(isnan([DCM.H{2}(:); DCM.H{1}(:)]))
            count=count+1
            nan_subs{count}=subjects{ss}
        end
    end
end

subjects_orig=subjects;
clearvars subjects

count=0;
for ss=1:length(subjects_orig)
    if ~any(find(contains(subjects_orig{ss},[failed_subs, nan_subs])))
        count=count+1;
        
        subjects{count}=subjects_orig{ss};
    end
end

BGCM=spm_dcm_load(BLbmr.RCM(:,maxind)); %maxind=6;

for f=1:length(BGCM)
    Bfiles{f}=BGCM{f}.name;
end

Bsubs=extractBetween(Bfiles, 'DCM_', '_')
Lsubs = subjects(ismember(subjects, Bsubs)); %remove any subs not in BL who are in AF (failed to converge)
save([anaL '/Lsubs'], 'Lsubs');

%% Plot model fits (line graphs)

if plot1 == 1
    
    dcms = cellstr(spm_select('FPList', [anaL filesep], '^*full*.mat'));
    dcms = spm_dcm_load(dcms);
    figure
    title('Model fit');
    
    for ss = 1:length(subjects)
        try
            DCM =spm_dcm_load([anaL filesep 'DCM_' subjects{ss} '_full.mat']);%(dcms{ss});
            DCM = DCM{1};
            reps=1;
            
            subplot(6,14,ss)
            
            for c = 1 %:5:size((DCM.H),2)
                plot(DCM.H{c}(:,1) + DCM.R{c}(:,1), 'color', [0.4 0.4 0.4], 'Linewidth', 1.5); hold on
                plot(DCM.H{c}(:,1), 'color', [0.3    0.8    0.6510], 'Linewidth', 1.5);
                
                obs1(ss,:)=DCM.H{c}(:,1) + DCM.R{c}(:,1);
                prd1(ss,:)=DCM.H{c}(:,1);
                
                cortemp=corrcoef(DCM.H{c}(:,1) + DCM.R{c}(:,1), DCM.H{c}(:,1));
                corL.s(ss,1)=subjects(ss);
                corL.d(ss,1)=cortemp(2);
                
                ylim([-7 4]); yticks([-6 4]);
                xlim([0 125]); xticks([25 75]); %xticklabels({'0', '200'});
                xlabel(DCM.name(end-9:end-5))
                set(gcf, 'color', 'w');
                set(gcf, 'Position', [100 + 400*(reps-1) 100 400 800]);
            end
            box off
            %maxP=max(DCM.H{c}(:,1))+1;
            %title(['r=' num2str(round(cor1(ss,1),2))], 'FontWeight', 'normal', 'Position', [75 maxP 0], 'FontSize', 5)

        catch; continue
        end
    end
    
    legend({'DEV(pred)', 'DEV(obs)' , 'REP5(pred)', 'REP5(obs)'});
    %exportgraphics(gcf, [anaL 'Lmodelfits.png'], 'Resolution', 720) %only works in 2020

    %average plot
    subplot(1,6,1); plot(nanmean(obs1,1), 'color', [0.4 0.4 0.4], 'Linewidth', 1.5); hold on; plot(nanmean(prd1,1), 'color', [0.3    0.8    0.6510], 'Linewidth', 1.5);
        
    %ylim([-5 5]);
    xlim([0 125]); xticks([25 75]); %xticklabels({'0', '200'});
    xlabel('longitudinal group average'); box off
    set(gcf, 'color', 'w');
    set(gcf, 'Position', [100 + 400*(reps-1) 100 400 800]);
end




%% get free energy of nested models for bayesian model comparison ...

if bmc == 1
    
    [DCMa] = dcm_cmc_gen_model_space_v1();
        
    % BASELINE: For each subject, generate reduced submodels to run BMR over
    %--------------------------------------------------------------------------
    
    ss = 1; % index of BL and AF of the 'P' variable (BLsubs + AFsubs)
    
    for s = 1:length(Lsubs)
        
        count = 1;
        
        %Add sub-models with exponential and phasic repetition effect
        %--------------------------------------------------------------------------
        for m = 1:(length(DCMa)-1)
            clear DCM
            DCM = dcm_cmc_gen_v1(DCMa{m}, raw, anaB, Lsubs{s}, rawlast);      %mmn_bmr_gen_jl(spatial, DCMa{l}, Fbmr);
            P{ss,count} = DCM;
            count = count + 1;           
        end
        
        % Add sub-models with exponential repetition only
        %--------------------------------------------------------------------------
        
        for m = 1:(length(DCMa)-1)
            clear DCM
            DCM = dcm_cmc_gen_expo_v1(DCMa{m}, raw, anaL, Lsubs{s}, rawlast);      %mmn_bmr_gen_jl(spatial, DCMa{l}, Fbmr);
            P{ss,count} = DCM;
            count = count + 1;            
        end
                
        % Add the full model that has been inverted to the end
        %--------------------------------------------------------------------------
        P{ss,count} = spm_dcm_load([anaL filesep 'DCM_' Lsubs{s} '_full.mat']);
        P{ss,count} = P{ss,count}{1};
        
        ss = ss+1;
                        
    end
    
    % FOLLOW UP
    %--------------------------------------------------------------------------
    %Generate reduced submodels from the full model for each subject
    
    for s = 1:length(Lsubs)
        
        count = 1;
        
        %Add sub-models with exponential and phasic repetition effect
        %--------------------------------------------------------------------------
        for m = 1:(length(DCMa)-1)
            clear DCM
            DCM = dcm_cmc_gen_v1(DCMa{m}, raw, anaL, Lsubs{s}, rawlast);      %mmn_bmr_gen_jl(spatial, DCMa{l}, Fbmr);
            P{ss,count} = DCM;
            count = count + 1;
        end
        
        % exponential decay only
        for m = 1:(length(DCMa)-1)
            clear DCM
            DCM = dcm_cmc_gen_expo_v1(DCMa{m}, raw, anaL, Lsubs{s}, rawlast);      %mmn_bmr_gen_jl(spatial, DCMa{l}, Fbmr);
            P{ss,count} = DCM;
            count = count + 1;            
        end
                
        % Add the full model that has been inverted to the end
        %--------------------------------------------------------------------------
        P{ss,count} = spm_dcm_load([anaL filesep 'DCM_' Lsubs{s} '_full.mat']);
        P{ss,count} = P{ss,count}{1};
        
        ss = ss+1;
    end
    
    %% Run BMR
    
    % MANUALLY DELETE FAILED dcms - 7/P1010, P1021, P1055, C1011 P1020 (some completely
    %failed, some have NANs )
    
    [RCM, BMC, BMA] = spm_dcm_bmr(P);    
    save([anaL filesep 'RCM_BMC_BMA_BLAFonlysubs'], 'RCM', 'BMC', 'BMA'); %save([anahe filesep 'BMC_BMA_2'], 'BMC', 'BMA')
    
end



%% Plot results of Bayesian model comparison FREE ENERGY
%%==========================================================================


if plot2 ==1
    
    AFbmr=load([anaL filesep 'RCM_BMC_BMA_BLAFonlysubs.mat']); 
    
    %% Plot Bayesian model comparison results
    %%=========================================================  
    clearvars F Fm dF Fms Fsort dF maxind
    figure
    subplot(2,1,1)
    F = [AFbmr.BMC.F]';
    Fma = mean(F,1);
    %Fm = Fma(1:24)
    Fm=Fma(1:12)
    Fm = Fm - min(Fm);
    Fms(:,1)=Fm(1:6)';
    Fms(:,2)=Fm(7:12)';
    
    Fsort = sort(Fm);
    dF    = Fsort(end) - Fsort(end-1);
    maxind= find(Fm==max(Fm)); nextind=find(Fm==Fsort(end-1)); %NB see bar(Fm) for correct model of max and next max - EG is on LHS
    
    b = bar(Fms);
    %b(1).FaceColor= [.2 .6 .6]; b(2).FaceColor = [0.4 0.4 0.4];
    b(1).FaceColor= [1 1 1]; b(2).FaceColor = [0.4 0.4 0.4];

    
    yticks([0 60000])
    ylabel('Free energy');
    legend('EG','E', 'Orientation', 'horizontal');legend('boxoff');
    %legend('EG','E'); %legend('Location', 'bestoutside')
    box off
    
    Pm = softmax(Fm');
    Pms(:,1)=Pm(1:6)';
    Pms(:,2)=Pm(7:12)';
    subplot(2,1,2)
    b = bar(Pms)
%     b(1).FaceColor= [.2 .6 .6]; b(2).FaceColor = [0.4 0.4 0.4];
    b(1).FaceColor= [1 1 1]; b(2).FaceColor = [0.4 0.4 0.4];
        
    %     title(['Posterior model probabilities']);
    xlabel(strcat('winningmodel is :', extractAfter(AFbmr.RCM{1, maxind}.name,[subjects{1} '_'])),'Interpreter', 'none'); %Models: No input (1:6) vs input(7:12) to parietal | No lateral (1:3) vs lateral (4:6) | no-IPC-IFG/ IPC > IFG/ IFG> IPC ' );
    ylabel('poseterior model probability'); yticks([0 1])
    box off
    set(gcf, 'color', 'white')
    
    exportgraphics(gcf, [anaF, '/BMC_PvC.png'], 'Resolution', '720');
end