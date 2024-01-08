%% sets xticklabels after running spm_dcm_peb_fig_jl
%   PEBx - Name of PEB / RMA used in figure
%   DCM /GCM used for reference (the one used in your spm_dcm_peb (if using
%   a group DCM (or GCM) specify a signle DCM i.e. GCM{1} 

function [xticknames] = xticklabels_jl_peb(PEBx, GCM)

for x=1:length(PEBx.Pnames)
    xticknames{1, x}=pname_to_string(PEBx.Pnames{x}, GCM.Sname);
    xticknames{2, x}=strcat(extractBetween(xticknames{1,x}, 'matrix  ', 'to'), 'to', xticknames{1,x}(end-3:end));
    xticknames{3, x}=strcat(extractBefore(xticknames{1,x}, '-'), ':  ', extractBetween(xticknames{1,x}, 'matrix  ', 'to'), 'to', xticknames{1,x}(end-3:end));
end

%xticklabels([xticknames{2,:}])
end