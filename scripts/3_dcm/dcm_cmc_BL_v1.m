%% DCM using the cmc model to test baseline differences (patients vs controls)
%  1. Inverts the full DCM for each individual and check fits
%  2. Runs a bayesian model reduction to compare models to identify where
%  parietal node is connected to other regions

%% Set up environment
clearvars
E = cmc_environment;

%% Set up variables
scr = E.scr;
load([scr '/BLsubs.mat']); 
subjects = BLsubs;
anaB = E.anaB; %character array to DCM files location
raw = E.raw;  % character aray to base file directory (e.g. raw='/juliette/cmc/meg_data')
rawlast = '/BL/mmn/ffmraeMaffffdtsss.mat'; % file that has had all NTAD preprocessing steps applied

%% switches
invert = 0;
plot_line = 1;
bmc = 0;
plot2 = 0;

%% Set up  model space

% %create models from my model space with different A and C matrices
[DCMa] = dcm_cmc_gen_model_space_v1();

% %specify rest of parameters for the full DCM model
for ss=1:length(subjects)
    dcms{ss} = dcm_cmc_gen_v1(DCMa{end}, raw, anaB, subjects{ss}, rawlast);
end

dcms=dcms';
GCM=spm_dcm_load(dcms);

%% Invert first-level DCMs
% % subject that do not converge go into the 'failed' variable

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
    save([anaB '/failedsubs.mat'],failed_subs)
end



%% Remove failed or bad fits
% From the variable 'failed_subs' created during the inversion and from
% seeing which have NANs in the H matrix we remove those failed subjects

if ~exist('failed_subs', 'var')
    load([anaB '/failedsubs.mat'])
end

count=0;
for ss = 1:length(subjects)
    try
        DCM =spm_dcm_load([anaB filesep 'DCM_' subjects{ss} '_full.mat']);
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

%% Plot model fits (line graphs)

if plot_line == 1
    
    dcms = cellstr(spm_select('FPList', [anaB filesep], '^*full*.mat'));
    dcms = spm_dcm_load(dcms);
    figure
    title('Model fit');
    
    for ss = 1:length(subjects)
        try
            DCM =spm_dcm_load([anaB filesep 'DCM_' subjects{ss} '_full.mat']);
            DCM = DCM{1};
            reps=1;
            
            subplot(6,14,ss)
            
            for c = 1
                plot(DCM.H{c}(:,1) + DCM.R{c}(:,1), 'color', [0.4 0.4 0.4], 'Linewidth', 1.5); hold on
                plot(DCM.H{c}(:,1), 'color', [0.3    0.8    0.6510], 'Linewidth', 1.5);
                
                obs1(ss,:)=DCM.H{c}(:,1) + DCM.R{c}(:,1);
                prd1(ss,:)=DCM.H{c}(:,1);
                
                cortemp=corrcoef(DCM.H{c}(:,1) + DCM.R{c}(:,1), DCM.H{c}(:,1));
                corL.s(ss,1)=subjects(ss);
                corL.d(ss,1)=cortemp(2);
                
                ylim([-5 5]);
                xlim([0 125]); xticks([25 75]); %xticklabels({'0', '200'});
                xlabel(DCM.name(end-9:end-5))
                set(gcf, 'color', 'w');
                set(gcf, 'Position', [100 + 400*(reps-1) 100 400 800]);
            end
            box off
            
        catch; continue
        end
    end
    legend({'DEV(pred)', 'DEV(obs)'});
    
    %exportgraphics(gcf, [anaF 'BLfits.png'], 'Resolution', 720) %only works in 2020
    
    %% group average plot
    figure
    subplot(1,6,1); plot(nanmean(obs1,1), 'color', [0.4 0.4 0.4], 'Linewidth', 1.5); hold on; plot(nanmean(prd1,1), 'color', [0.3    0.8    0.6510], 'Linewidth', 1.5);
    cortemp=corrcoef(nanmean(prd1,1), nanmean(obs1,1));
    corav=cortemp(2);
    title(['r=' num2str(round(corav,2))], 'FontWeight', 'normal', 'Position', [75 0.8, 0], 'FontSize', 6)
    
    ylim([-5 5]);
    xlim([0 125]); xticks([25 75]); %xticklabels({'0', '200'});
    xlabel('baseline group average'); box off
    set(gcf, 'color', 'w');
    set(gcf, 'Position', [100 + 400*(reps-1) 100 400 800]);
end

%% get free energy of reduced models for bayesian model comparison

if bmc == 1
    
    [DCMa] = dcm_cmc_gen_model_space_v1();
    
    %  Generate reduced submodels from the full model for each subject
    
    for ss = 1:length(subjects)
        
        count = 1;
        
        % Add sub-models with exponential and phasic repetition effect
        
        for m = 1:(length(DCMa)-1)
            clear DCM
            DCM = dcm_cmc_gen_v1(DCMa{m}, raw, anaB, subjects{ss}, rawlast);      
            P{ss,count} = DCM;
            count = count + 1;
            
        end
        
        % Add sub-models with exponential repetition only (so B{2}=
        % zero)
        
        for m = 1:(length(DCMa)-1)
            clear DCM
            DCM = dcm_cmc_gen_expo_v1(DCMa{m}, raw, anaB, subjects{ss}, rawlast);    
            P{ss,count} = DCM;
            count = count + 1;
        end
        
        % Add the full model that has been inverted to the end
        P{ss,count} = spm_dcm_load([anaB filesep 'DCM_' subjects{ss} '_full.mat']);
        P{ss,count} = P{ss,count}{1};
        
    end
    
    
    %  Run BMR
    [RCM, BMC, BMA] = spm_dcm_bmr(P);
    
    save([anaB filesep 'RCM_BMC_BMA_may22'], 'RCM', 'BMC', 'BMA'); 
end

if plot2 ==1
    
    % Plot Bayesian model comparison results
    
    if ~exist('BMC', 'var')
        load([anaB filesep 'RCM_BMC_BMA_may22']);
    end
    
    clearvars F Fm dF Fms Fsort dF maxind
    figure
    subplot(2,1,1)
    F = [BMC.F]';
    Fma = mean(F,1);
    Fm=Fma(1:12)
    Fm = Fm - min(Fm);
    Fms(:,1)=Fm(1:6)';
    Fms(:,2)=Fm(7:12)';
    
    Fsort = sort(Fm);
    dF    = Fsort(end) - Fsort(end-1);
    maxind= find(Fm==max(Fm)); nextind=find(Fm==Fsort(end-1)); %NB see bar(Fm) for correct model of max and next max - EG is on LHS
    
    b = bar(Fms);
    b(1).FaceColor= [.2 .6 .6]; b(2).FaceColor = [0.4 0.4 0.4];
    b(1).FaceColor= [1 1 1]; b(2).FaceColor = [0.4 0.4 0.4];
    
    yticks([0 60000])
    ylabel('Free energy');
    legend('EG','E', 'Orientation', 'horizontal','Location', 'northoutside'); legend('boxoff');
    
    box off
    
    
    Pm = softmax(Fm');
    Pms(:,1)=Pm(1:6)';
    Pms(:,2)=Pm(7:12)';
    subplot(2,1,2)
    b = bar(Pms)
    b(1).FaceColor= [.2 .6 .6]; b(2).FaceColor = [0.4 0.4 0.4];
    b(1).FaceColor= [1 1 1]; b(2).FaceColor = [0.4 0.4 0.4];
    
    
    %title(['Posterior model probabilities']);
    xlabel(strcat('winningmodel is :', DCMa{maxind}.name),'Interpreter', 'none'); %Models: No input (1:6) vs input(7:12) to parietal | No lateral (1:3) vs lateral (4:6) | no-IPC-IFG/ IPC > IFG/ IFG> IPC ' );
    ylabel('poseterior model probability'); yticks([0 1]);
    box off
    set(gcf, 'color', 'white')
    
   % exportgraphics(gcf, [scr filesep 'BMC_PvC.png'], 'Resolution', '720');
end