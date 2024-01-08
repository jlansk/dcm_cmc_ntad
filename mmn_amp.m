function [mmn_pa, mmn] = mmn_amp(subjects, ana_dir, session, task, lfile, erf)
% Calculate the mmn amplitude
%   subjects: 1 x n cell of subject IDs
%   ana_dir: string to file location of the meg data folder
%   session: recording session ('BL' / 'AF' / 'TW')
%   task: MMN task ('MMN')
%   lfile: name of file 
%   erf: whether to generate data for MMN or ERF (for MMN = 0, ERF = 1)

%% Load the data

% % load data
for ss=1:length(subjects)
    ptss=load(strcat(ana_dir, subjects{ss}, filesep, session, filesep, task, filesep, lfile));
    mpind=find(contains({ptss.D.channels.type}, 'MEGCOMB'));
    ptsss(ss,:,:)=nanmean(ptss.D.data(mpind,:,:),1);
end
data=double(ptsss);

%% Calculate the amplitude of the mismatch negativity waveform

for s=1:size(data,1) %for each sub
    for c=2:size(data,3)
        mmn(s, :, c-1)=squeeze(data(s,:,c)-data(s,:,1)); %standard minus oddball
    end
end

if erf == 1
    mmn = data;
end

%% find peak MMN waveform amplitude

for s=1:size(mmn,1) %for each subject
    
    begin2=140; beg2=(begin2+100)/2; %significant diff between rep6 and dev for whole group (bonferonni corrected)
    finish2=160 ; fin2=(finish2+100)/2;
    
    for c=1:size(mmn,3)
        data2=(squeeze(mmn(s,beg2:fin2,c)));
        if any(data2(1))
            mmn_pa((c), s) = nanmean(data2);
        end
    end
end

mmn_pa=mmn_pa';
end
