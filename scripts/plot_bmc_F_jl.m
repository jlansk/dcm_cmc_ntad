function plot_bmc_F_jl(F, labels, title_extra)
% Generates a plot with differential free energies across models and
% posterior model probabilities.
% F should be double with free energies in 1 x (Number of models) double
% (e.g. F(1)=PEB(1).F; F(2)=PEB(2).F). labels should be cell with string of
% model names for xaxis. title_extra can be additional info for title if
% needed and should be character array e.g. title_extra=('Group BMC');

% Prepare input
% -------------------------------------------------------------------------
if nargin < 3 || isempty(title_extra)
    title_extra=0;

end

figure
subplot(2,1,1)
b=bar(F-min(F));
b.FaceColor=[0.4 0.4 0.4];

% Labels and Titles
Fsort = sort(F, 'descend');
df = Fsort(2) - Fsort(1);
ylabel('Relative Free Energy');

if title_extra
    title([title_extra ', \DeltaF = ' num2str(df)], 'Fontsize', 15);
else
    title(['\DeltaF = ' num2str(df)], 'Fontsize', 15);
end
box off

% Plot settings
xticks([1:1:length(F)]);
set(gca, 'XTickLabels', labels, 'FontSize', 12);
xtickangle(45)

% Plot posterior model probability
%--------------------------------------------------------------------------
subplot(2,1,2)
title('Posterior Probability');
b=bar(spm_softmax(F'));
b.FaceColor=[0.4 0.4 0.4];

ylabel('Posterior Probability');

set(gcf, 'color', 'w');
set(gcf, 'Position', [100 500 900 400]);
box off
xticks([1:1:length(F)]);
set(gca, 'XTickLabels', labels, 'FontSize', 12);

set(gcf, 'units','normalized','outerposition',[0 0 1 1])
       xtickangle(45) 
end
