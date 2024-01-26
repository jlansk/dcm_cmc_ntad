%% Sensor level analysis - Juliette Lanskey 2023
%This script runs the sensor space t-tests and ANCOVAs for the baseline
%(controls vs AD/MCI) and longitudinal (AD/MCI baseline vs follow up)
%comparisons and generates mismatch negativity waveform plots.

%v2 - longitudinal ANOVA adjusted to make session within-subject rather than between (as it is
%a longitudinal ANOVA)


% set up variables
clearvars
E = cmc_environment;
raw = E.raw; %directory for MEG files
scr = E.scr; % directory of scripts
ana_dir = raw; % where the pre-processed MEG files are

load([scr filesep 'BLsubs.mat']); % 1x59 cell with baseline IDs
load([scr filesep 'AFsubs.mat']); % 1x33 cell with follow-up IDs

task='mmn';
session1 = 'BL';
session2 = 'AF';
lfile= 'bCPffmraeMaffffdtsss.mat'; % file with preprocessing and combined grads

%% load baseline MMN amplitude
[mmn_pa_BL, mmnBL] = mmn_amp(BLsubs, ana_dir, 'BL', task, lfile, 0);

%% load follow-up MMN amplitude
[mmn_pa_AF, mmnAF] =mmn_amp(AFsubs, ana_dir, 'AF', task, lfile, 0); %12 is dev-stn, not needed

%% statistics

%Baseline Group comparison
cons = find(contains(BLsubs, 'C'));
pats = find(contains(BLsubs, 'P'));

[h,p,ci,stats]=ttest2(mmn_pa_BL(cons,1),mmn_pa_BL(pats,1), 'tail','left');
cnanmean=nanmean(mmn_pa_BL(cons,1)); pnanmean=nanmean(mmn_pa_BL(pats,1));
BLstatistics(1,1)=stats.tstat; BLstatistics(2,1)=p; BLstatistics(3,1)=cnanmean; BLstatistics(4,1)=pnanmean;

%Longitudinal comparison
pats = find(contains(BLsubs, 'P'));

Lsubs = AFsubs(ismember(AFsubs, BLsubs));
LBsubs = find(contains(BLsubs, Lsubs));
LAsubs = find(contains(AFsubs, Lsubs));


[h,p,ci,stats]=ttest(mmn_pa_BL(LBsubs,1),mmn_pa_AF(LAsubs,1), 'tail','left');
bnanmean=nanmean(mmn_pa_BL(LBsubs,1)); fnanmean=nanmean(mmn_pa_AF(LAsubs,1));

AFstatistics(1,1)=stats.tstat; AFstatistics(2,1)=p; AFstatistics(3,1)=bnanmean; AFstatistics(4,1)=fnanmean;


%% ANOVAs
mmn_pa_BL_ERF = mmn_amp(BLsubs, ana_dir, 'BL', task, lfile, 1);
group_pa_BL(cons,1)=0; group_pa_BL(pats,1)=1;

mmn_pa_AF_ERF=mmn_amp(AFsubs, ana_dir, 'AF', task, lfile, 1);

%% ANOVA with group, tone and age

load([scr filesep 'BLage.mat']);

Bsubs=BLsubs
for ss=1:length(Bsubs)
    idx=find(contains(BLage(:,1),Bsubs(ss)));
    Xage(ss,1) = BLage(idx,1);
    Xage(ss,2) = BLage(idx,2); % should be the value corresponding to the peb_subs
end

%With dev:rep5
group_pa_BL(cons,1)=0; group_pa_BL(pats,1)=1;
t = table(mmn_pa_BL_ERF(:,1), mmn_pa_BL_ERF(:,2), mmn_pa_BL_ERF(:,3), mmn_pa_BL_ERF(:,4), mmn_pa_BL_ERF(:,5), mmn_pa_BL_ERF(:,6), group_pa_BL(:,1), zscore([Xage{:,2}]'), 'VariableNames', {'rep0', 'rep1', 'rep2', 'rep3', 'rep4', 'rep5', 'group', 'age'});
trial = table([1 2 3 4 5 6]', 'VariableNames', {'trial'});
rm = fitrm(t, 'rep0-rep5~group+age', 'WithinDesign', trial);
ranovatbl=ranova(rm);

%% Longitudinal
mmn_pa_BF = [mmn_pa_BL_ERF(LBsubs,1:6), mmn_pa_AF_ERF(LAsubs,1:6)];

t2 = table(mmn_pa_BF(:,1), mmn_pa_BF(:,2), mmn_pa_BF(:,3), mmn_pa_BF(:,4), mmn_pa_BF(:,5), mmn_pa_BF(:,6), mmn_pa_BF(:,7), mmn_pa_BF(:,8), mmn_pa_BF(:,9), mmn_pa_BF(:,10), mmn_pa_BF(:,11), mmn_pa_BF(:,12),'VariableNames', {'Brep0', 'Brep1', 'Brep2', 'Brep3', 'Brep4', 'Brep5', 'Lrep0', 'Lrep1', 'Lrep2', 'Lrep3', 'Lrep4', 'Lrep5'});
WithinStruct = table([1 1 1 1 1 1 2 2 2 2 2 2]', [1 2 3 4 5 6 1 2 3 4 5 6]', 'VariableNames', {'Session', 'Trial'});
WithinStruct.Treatment = categorical(WithinStruct.Session);
WithinStruct.Trial = categorical(WithinStruct.Trial);

rm = fitrm(t2, 'Brep0,Brep1,Brep2,Brep3,Brep4,Brep5,Lrep0,Lrep1,Lrep2,Lrep3,Lrep4,Lrep5~1', 'WithinDesign', WithinStruct);
ranovatable = ranova(rm, 'WithinModel', 'Session*Trial');

%% bar charts - June 22
% Baseline control vs patients average ERF between 140-160ms for all tones
CP_pa(:,1) = nanmean(mmn_pa_BL_ERF(cons,1:11),1);
CP_pa(:,2) = nanmean(mmn_pa_BL_ERF(pats,1:11),1);
figure(1)
b=bar(CP_pa);
b(1).FaceColor=[1 0.5 0.16]; %[0.9 0.7 0.3];%[0.55 0.55 1];
b(2).FaceColor=[0 0.32 1]; %[0.9 0.7 0.3];
box off; set(gcf, 'color', 'w'); %ylim([0 1.5]);

%BL vs AF bar plots for erfs
figure(2)
BF_pa(:,1) = nanmean(mmn_pa_BL_ERF(LBsubs,1:11),1);
BF_pa(:,2) = nanmean(mmn_pa_AF_ERF(LAsubs,1:11),1);
b=bar(BF_pa);
b(1).FaceColor=[0 0.32 1]; %[0.9 0.7 0.3];
b(2).FaceColor=[0.05 0.8 0.85]; %[0.8 0.3 0.2];
box off; set(gcf, 'color', 'w'); %ylim([0 1.5]);

%% plots

% MMN waveforms
%if mmn1 ->

BL_av = squeeze(nanmean(mmnBL(LBsubs,:,:),1))
BL_se= squeeze(nanstd(mmnBL(LBsubs,:,:),1))./sqrt(size(LBsubs,2));

AF_av = squeeze(nanmean(mmnAF(LAsubs,:,:),1))
AF_se = squeeze(nanstd(mmnAF(LAsubs,:,:),1))./sqrt(size(LAsubs,2))


con_av = squeeze(nanmean(mmnBL(cons,:,:),1))
con_se=squeeze(nanstd(mmnBL(cons,:,:),1))./sqrt(size(cons,2))%standard error

pat_av = squeeze(nanmean(mmnBL(pats,:,:),1))
pat_se = squeeze(nanstd(mmnBL(pats,:,:),1))./sqrt(size(pats,2))

%% repetition plots seperately and saved for aruk east

titlestr={'MMN_r1r0'};
choice = menu('Save figures?', 'Yes', 'No');


%% BL AD VS HC
close all

bh(1)=plot(con_av(:,1), 'Color', [1 0.42 0.16]); hold on; bh(2)=plot(pat_av(:,1), 'Color', [0 0.32 1]); hold on;
set(bh(1), 'linewidth', 3); set(bh(2), 'linewidth', 3);
boundedline([1:size(con_av)],con_av(:,1),con_se(:,1), 'transparency', 0.4, 'cmap', [1 0.42 0.16]); boundedline([1:size(pat_av)],pat_av(:,1),pat_se(:,1), 'transparency', 0.4, 'cmap', [0 0.1 0.75]); %see uisetcolor for colour options
xlim([0 250]); xticks([0 50 100 150 200 250]); xticklabels({'-100','0', '100', '200', '300', '400'}); xlabel('Time (ms)');
xlabel('Time (ms)'); ylabel('Mismatch response (fT/m)');
begin2=140; beg2=(begin2+100)/2;
finish2=160 ; fin2=(finish2+100)/2;
ylim([-1 0.1])
% patch([beg2,beg2, fin2, fin2, beg2], [0.1, -1.0, -1.0, 0.1, 0.1], 'black', 'EdgeColor', 'none', 'FaceColor', 'black', 'FaceAlpha', 0.1);
box off
set(gcf, 'color', 'w')

legend('controls', 'patients'); legend('Location', 'best')

if choice ==2 | choice ==0
    return;
elseif choice ==1
    try
        exportgraphics(gcf, [scr filesep 'CP_' titlestr{c} '.png'], 'Resolution', 720) %only works in 2020
    catch
        saveas(gca, [scr filesep 'CP_' titlestr{c} '.fig'])
    end
end

%% longitudinal plot
close all
lh(1)=plot(BL_av(:, 1), 'Color', [0 0 0.6]); hold on; lh(2)=plot(AF_av(:, 1), 'Color', [0.05 0.6 0.65]); hold on;
set(lh(1), 'linewidth', 3); set(lh(2), 'linewidth', 3);
boundedline([1:size(BL_av)],BL_av(:,1),BL_se(:,1), 'transparency', 0.333, 'cmap', [0 0.1 0.7]); boundedline([1:size(AF_av)],AF_av(:,1),AF_se(:,1), 'transparency', 0.333, 'cmap',  [0.05 0.6 0.75]); %[0.05 0.8 0.85]see uisetcolor for colour options
xlim([0 250]); xticks([0 50 100 150 200 250]); xticklabels({'-100','0', '100', '200', '300', '400'}); xlabel('Time (ms)');
xlabel('Time (ms)'); ylabel('Mismatch response (fT/m)');
%  patch([beg2,beg2, fin2, fin2, beg2], [0.1, -1.0, -1.0, 0.1, 0.1], 'black', 'EdgeColor', 'none', 'FaceColor', 'black', 'FaceAlpha', 0.1);
begin2=140; beg2=(begin2+100)/2; 
finish2=160 ; fin2=(finish2+100)/2;
ylim([-1 0.1]); box off; set(gcf, 'color', 'w');
legend('baseline', 'follow-up'); legend('Location', 'best')

if choice ==2 | choice ==0
    return;
elseif choice ==1
    try
        exportgraphics(gcf, [scr filesep '/BF_' titlestr{1} '.png'], 'Resolution', 720) %only works in 2020
    catch
        saveas(gca, [scr filesep '/BF_' titlestr{1} '.fig'])
    end
    
end
    
