%% Test retest
% Calculates the absolute intraclass correlation for the baseline and
% two-week MEG scans

% set up environment and variables
clear
E = cmc_environment; 
scr = E.scr;
raw = E.raw;

load([scr filesep 'TWsubs.mat']) % 1x14 cell with test-retest IDs
subjects=TWsubs;
task='mmn';
lfile= 'bCPffmraeMaffffdtsss.mat';
ana_dir = raw;
basefile = 'tsss';
condnames = {'dev','rep1', 'rep2','rep3','rep4', 'rep5', 'rep6', 'rep7', 'rep8', 'rep9', 'rep10'};


%% Calculate the MMN waveform amplitude
mmn_pa_BL = mmn_amp(subjects, ana_dir, 'BL', task, lfile);
mmn_pa_TW = mmn_amp(subjects, ana_dir, 'TW', task, lfile);

%% -------------------------------------------------
%%  ICCs - run above first for BL and TW

data = mmn_pa_BL(any(mmn_pa_BL,2),:); data = data(:,any(data,1)); data = data(any(data,2),:);

for s=1:size(data,2) %for each condition
    
    for r=1:size(data,1) %for each subject
        
        if ~(data(r,s)==0)
            mmn_pa_icc(size(data,1)*s+r,1)=data(r,s);%should be size conditions*pats + pats
        end
    end
end

data = mmn_pa_TW(any(mmn_pa_TW,2),:); data = data(:,any(data,1)); data = data(any(data,2),:);

for s=1:size(data,2)
    
    for r=1:size(data,1)
        
        if ~(data(r,s)==0)
            mmn_pa_icc(size(data,1)*s+r,2)=data(r,s);
        end
    end
end

mmn_pa_icc(1:(size(data,1)),:)=[]

allcondnames = {'all_mmn1', 'all_mmn2', 'all_mmn3', 'all_mmn4','all_mmn5','all_mmn6','all_mmn7','all_mmn8','all_mmn9','all_mmn10'};

for r=1:(size(data,2))

    data1=mmn_pa_icc(((r*length(subjects)-(length(subjects)-1)):r*length(subjects)),1:2);
    data1 = data1(all(data1,2),:);  %data1 is the BL (collumn 1) and TW (collumn 2) values for allcondame(r) - it deletes any subs who had a 0 at one or both sessions
    
    [R, LB, UB, F, df1, df2, p] = ICC(data1,'A-k')
    
    icc_r_pa(2,1) = {'mmn_pa_r'};
    icc_r_pa(3,1) = {'mmn_pa_LB'};
    icc_r_pa(4,1) = {'UB'};
    icc_r_pa(5,1) = {'F'};
    icc_r_pa(6,1) = {'df1'};
    icc_r_pa(7,1) = {'df2'};
    icc_r_pa(8,1) = {'p'};
    
    icc_r_pa(1,r+1) = {allcondnames{r}};
    icc_r_pa(2,r+1) = {R};
    icc_r_pa(3,r+1) = {LB};
    icc_r_pa(4,r+1) = {UB};
    icc_r_pa(5,r+1) = {F};
    icc_r_pa(6,r+1) = {df1};
    icc_r_pa(7,r+1) = {df2};
    icc_r_pa(8,r+1) = {p};
    
end

%% scatter plot
r=1;

data1=mmn_pa_icc(((r*length(subjects)-(length(subjects)-1)):r*length(subjects)),1:2);
data1 = data1(all(data1,2),:);  %data1 is the BL (collumn 1) and TW (collumn 2) values for allcondame(r) - it deletes any subs who had a 0 at one or both sessions

[R, LB, UB, F, df1, df2, p] = ICC(data1,'A-k');

x=data1(:,1);
y=data1(:,2);

scatter(x,y, 'filled', 'MarkerFaceColor', [0.5 0.5 0.5]); xlabel('Baseline MMN amplitude (fT/m)'); ylabel('Two-week retest MMN amplitude (fT/m)');title(['Test retest absolute intraclass correlation (R=', num2str(R), 'p=', num2str(p), ')']);%<.001)');
xlim([-1.5 0.2]); ylim([-1.5 0.2]);
axis square
c=polyfit(x,y,1); y_est=polyval(c,x); set(gcf, 'color', 'w'); box off;
hold on;


%exportgraphics(gcf, [scr '/figs/ICC.png'], 'Resolution', 720) %only works in 2020