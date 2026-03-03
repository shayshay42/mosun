clc
clear all;

% Run either of these scenarios
cases = [10, 11, 12, 13]; % Preclincal Study 1 
cases = [20, 21, 22]; % Preclincal Study 2
cases = [30]; % Blin ALL Study


cases = [10, 11, 12, 13, 20, 21, 22, 30];

for case_no_dbl = cases
    close all;
    case_no_int = floor(case_no_dbl)
    clearvars -except case_no_int case_no_dbl color_str h_plot cases G_cnt;
   
    s = sbioloadproject('TDBr26_6_paper.sbproj');
    c = struct2cell(s);
    clear s;
    model = c{1};
    
    switch case_no_int
        case {10, 11, 12, 13}
            variant_id = [5; 9; 11; 13; 14; 24];
        case {20, 21, 22}
            variant_id = [5; 9; 11; 13; 14; 15; 24];
        case {30}
            variant_id = [5; 9; 11; 13; 14; 20; 22; 24];
    end
    
    for ii = 1 : length(variant_id)
        model.variant(variant_id(ii)).Active = true;
    end
    
    switch case_no_int
        case {10, 11, 12, 13, 20, 21, 22}
            dose_id = [7, 8, 18, 19];
        case {30}
            dose_id = [8, 19, 20];   
    end
    for ii = 1 : length(dose_id)
        dose_names{ii} = model.dose(dose_id(ii)).Name;
    end
    
    
    switch case_no_int
        case {10, 11, 12, 13}
            outvec = {'Bpb_perml'; 'totTpb_perml'; 'Tafraction_pb'; 'BTtotRatio_tiss'; 'BTtotRatio_tiss2'; 'IL6combo'}; 
        case {20, 21, 22} 
            outvec = {'Bpb_perml'; 'totTpb_perml'; 'Tafraction_pb';};
        case {30} 
            outvec = {'Bpb_perml'; 'totTpb_perml'; 'Tafraction_pb'};
    end
    
    switch case_no_int
        case {10, 11, 12, 13}
            filename='dataset_1.xlsx';
        case {20, 21, 22}
            filename='dataset_2.xlsx';
        case {30}
            filename='dataset_blin.xlsx'; 
    end
    [data]=xlsreadstruct(filename);
    
    switch case_no_int
        case {10, 11, 12, 13}
            datavec = {'Bcells'; 'CD8Tcells'; 'fracCD69'; 'spl_B_CD8'; 'LN_combined_B_CD8'; 'IL6'};
        case {20, 21, 22}
            datavec = {'Bcells'; 'CD8Tcells'; 'fracCD69';}; 
        case {30}
            datavec = {'Bcells'; 'CD8Tcells'; 'fracCD69'};
    end
    
    switch case_no_int
        case {10, 11, 12, 13, 20, 21, 22}
            paramNames = [[], {'Bpbo_perml', 'Bpbref_perml', 'Trpbo_perml', 'Trpbref_perml', 'PKflag', 'VPid', 'end_time'}];  
        case {30}
            paramNames = [[], {'Bpbo_perml', 'Bpbref_perml', 'Trpbo_perml', 'Trpbref_perml', 'PKflag'}]; 
    end
    
    for ii = 1 : length(paramNames)
        parameters(ii) = sbioselect(model, 'Name', paramNames{ii});
    end
    
    Tdata_tmp = sort(unique(data.Time));
    SimTime = max(Tdata_tmp) + 0.1;
    outtime = sort(unique([(0:0.1:SimTime)';Tdata_tmp]));
    outtime = outtime(outtime >= 0);
    
    expmodel = export(model, parameters);
    accelerate(expmodel);
    expmodel.SimulationOptions.OutputTimes = outtime;
    expmodel.SimulationOptions.StopTime = SimTime;
    
    [h_plot(case_no_int+1, :), G_cnt, Obj] = objective_run(case_no_int, [], [], expmodel, data, datavec, outvec, dose_names, 'g');
end




