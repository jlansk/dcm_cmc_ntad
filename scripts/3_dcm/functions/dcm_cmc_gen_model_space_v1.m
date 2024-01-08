function [DCMa] = dcm_cmc_gen_model_space_v1()

% Generating the alternate model structures
%==========================================================================
%creates A and C matrices from model space for Nareas = 8. The connections
%and inputs involving IPC differ across models.

% first represent areas as an index to make specifying matrices easier

LA1=1; 
LSTG=2;
LIFG=3;
LIPC=4;

RA1=5;
RSTG=6;
RIFG=7;
RIPC=8;

Nareas=8;

for da=1:3 % create A matrices of first 3 models
    
    %these models have no intrinsic or input to parietal

        DCMa{da}.A{1} = zeros(Nareas, Nareas);   % forward connections
        DCMa{da}.A{1}(LSTG,LA1) = 1; %to STG from A1
        DCMa{da}.A{1}(LIFG,LSTG) = 1; %to IFG from STG
        DCMa{da}.A{1}(LIPC,LSTG) = 1; %to IPC from STG
        
        DCMa{da}.A{1}(RSTG,RA1) = 1;
        DCMa{da}.A{1}(RIFG,RSTG) = 1;
        DCMa{da}.A{1}(RIPC,RSTG) = 1;
   
        if da==2
            
            DCMa{da}.A{1}(LIPC, LIFG) = 1;
            DCMa{da}.A{1}(RIPC,RIFG) = 1;
        
        elseif da==3
            
               DCMa{da}.A{1}(LIFG,LIPC) = 1;
               DCMa{da}.A{1}(RIFG,RIPC) = 1;
        end
        
        DCMa{da}.A{2} = zeros(Nareas, Nareas);    % backward connections
        DCMa{da}.A{2}(LA1,LSTG) = 1; %to A1 from STG
        DCMa{da}.A{2}(LSTG, LIFG) = 1;
        DCMa{da}.A{2}(LSTG,LIPC) = 1;

        DCMa{da}.A{2}(RA1,RSTG) = 1;
        DCMa{da}.A{2}(RSTG,RIFG) = 1;
        DCMa{da}.A{2}(RSTG,RIPC) = 1;
        
        if da==2
            
            DCMa{da}.A{2}(LIFG,LIPC) = 1;
            DCMa{da}.A{2}(RIFG,RIPC) = 1;
        
        elseif da==3
            
               DCMa{da}.A{2}(LIPC, LIFG) = 1;
               DCMa{da}.A{2}(RIPC,RIFG) = 1;
        end

           
        DCMa{da}.A{3} = zeros(Nareas,Nareas);    % modulatory connections/A{3} binary constraints on extrinsic connections
        DCMa{da}.A{3}(LA1, LA1) = 1;
        DCMa{da}.A{3}(LSTG,LSTG) = 1;
        DCMa{da}.A{3}(LIFG, LIFG) = 1;
        DCMa{da}.A{3}(LIPC,LIPC) = 1;
        DCMa{da}.A{3}(RA1,RA1) = 1;
        DCMa{da}.A{3}(RSTG,RSTG) = 1;
        DCMa{da}.A{3}(RIFG,RIFG) = 1;
        DCMa{da}.A{3}(RIPC,RIPC) = 1;

        
end

for da=4:6 %duplicate generated A matrices 
    DCMa{da}=DCMa{da-3};
end

for da =1:3 %specify C matrices (first 6)
    DCMa{da}.C = [1; 0; 1; 0; 1; 0; 1; 0];
end

for da=4:6 %specify C matrices (last 6)
    DCMa{da}.C= [1; 0; 1; 1; 1; 0; 1; 1]; 
end



%full model
DCMa{7}=DCMa{6};

for a =1:length(DCMa{5}.A)
    DCMa{7}.A{a}=logical(DCMa{6}.A{a}+DCMa{5}.A{a});
end


%% Names

% No input to parietal

% % % No connections between parietal and frontal
DCMa{1}.name = 'no-hl_no-pi';

% % % parietal higher than frontal
DCMa{2}.name = 'p-f_no-pi';

% % % frontal higher than parietal
DCMa{3}.name = 'f-p_no-pi';

% input to parietal
DCMa{4}.name = 'no-hl_int_pi';
DCMa{5}.name = 'p-f_int_pi';
DCMa{6}.name = 'f-p_int_pi';

DCMa{7}.name = 'full'

    
