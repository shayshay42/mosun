set(0, 'DefaultFigureVisible', 'off');

repo_root = fileparts(fileparts(mfilename('fullpath')));
out_dir = fullfile(repo_root, 'generated', 'matlab_reference');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

cases = [10, 20, 30];
ids = [2001, 1001, 1];
meta_rows = {};
row = 1;

for k = 1:length(cases)
    [sim_key, case_no, id_val, schedule, outvec, paramNames, paramValues, selectedDoseNames, Tsim, Xsim] = run_case(cases(k), ids(k));

    out_csv = fullfile(out_dir, [sim_key, '.csv']);
    out_tbl = array2table([Tsim, Xsim], 'VariableNames', ['time', outvec]);
    writetable(out_tbl, out_csv);

    rel_csv = fullfile('generated', 'matlab_reference', [sim_key, '.csv']);
    meta_rows{row,1} = sim_key;
    meta_rows{row,2} = case_no;
    meta_rows{row,3} = id_val;
    meta_rows{row,4} = schedule;
    meta_rows{row,5} = jsonencode(outvec);
    meta_rows{row,6} = jsonencode(paramNames);
    meta_rows{row,7} = jsonencode(paramValues);
    meta_rows{row,8} = jsonencode(selectedDoseNames);
    meta_rows{row,9} = rel_csv;
    row = row + 1;

    fprintf('Saved %-16s | ntime=%4d | schedule=%s\n', sim_key, length(Tsim), schedule);
end

meta_tbl = cell2table(meta_rows, 'VariableNames', {
    'sim_key', 'case_no', 'id', 'schedule', 'outvec_json', ...
    'param_names_json', 'param_values_json', 'selected_doses_json', 'sim_csv_relpath'});
writetable(meta_tbl, fullfile(out_dir, 'manifest_three_cases.csv'));

fprintf('Wrote manifest_three_cases.csv\n');


function [sim_key, case_no, id_val, schedule, outvec, paramNames, param_values, selected_dose_names, Tsim, Xsim] = run_case(case_no, id_val)
    matlab_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'assets', 'Supp Matlab Code');
    cd(matlab_dir);

    s = sbioloadproject('TDBr26_6_paper.sbproj');
    c = struct2cell(s);
    model = c{1};

    switch floor(case_no)
        case {10, 11, 12, 13}
            variant_id = [5; 9; 11; 13; 14; 24];
            dose_id = [7, 8, 18, 19];
            filename = 'dataset_1.xlsx';
            paramNames = {'Bpbo_perml', 'Bpbref_perml', 'Trpbo_perml', 'Trpbref_perml', 'PKflag', 'VPid', 'end_time'};
            outvec = {'Bpb_perml', 'totTpb_perml', 'Tafraction_pb', 'BTtotRatio_tiss', 'BTtotRatio_tiss2', 'IL6combo'};
        case {20, 21, 22}
            variant_id = [5; 9; 11; 13; 14; 15; 24];
            dose_id = [7, 8, 18, 19];
            filename = 'dataset_2.xlsx';
            paramNames = {'Bpbo_perml', 'Bpbref_perml', 'Trpbo_perml', 'Trpbref_perml', 'PKflag', 'VPid', 'end_time'};
            outvec = {'Bpb_perml', 'totTpb_perml', 'Tafraction_pb'};
        case {30}
            variant_id = [5; 9; 11; 13; 14; 20; 22; 24];
            dose_id = [8, 19, 20];
            filename = 'dataset_blin.xlsx';
            paramNames = {'Bpbo_perml', 'Bpbref_perml', 'Trpbo_perml', 'Trpbref_perml', 'PKflag'};
            outvec = {'Bpb_perml', 'totTpb_perml', 'Tafraction_pb'};
        otherwise
            error('Unsupported case_no: %d', case_no);
    end

    for ii = 1:length(variant_id)
        model.variant(variant_id(ii)).Active = true;
    end

    dose_names = cell(1, length(dose_id));
    for ii = 1:length(dose_id)
        dose_names{ii} = model.dose(dose_id(ii)).Name;
    end

    for ii = 1:length(paramNames)
        parameters(ii) = sbioselect(model, 'Name', paramNames{ii}); %#ok<SAGROW>
    end

    data = xlsreadstruct(filename);
    pos = find(data.ID == id_val);
    Tdata = sort(unique(data.Time(pos)));

    SimTime = max(Tdata) + 0.1;
    if case_no == 30
        SimTime = max(Tdata);
    end
    outtime = sort(unique([(0:0.1:SimTime)'; Tdata]));
    outtime = outtime(outtime >= 0);

    em = export(model, parameters);
    em.SimulationOptions.OutputTimes = outtime;
    em.SimulationOptions.StopTime = SimTime;

    init_Bpbo = mean(data.Bcells(pos(data.Time(pos)==0)), 'omitnan') * 1000;
    init_Tpbo = mean(data.CD8Tcells(pos(data.Time(pos)==0)), 'omitnan') * 1000;

    if case_no ~= 30
        param_values = [init_Bpbo, init_Bpbo, init_Tpbo, init_Tpbo, 0, id_val, SimTime];
    else
        param_values = [init_Bpbo, init_Bpbo, init_Tpbo, init_Tpbo, 1];
    end

    sch = unique(data.schedule(pos));
    if iscell(sch)
        schedule = char(sch{1});
    else
        schedule = char(sch);
    end

    if strcmp(schedule, '4qw')
        if any(case_no == [10, 11, 12, 13])
            if ismember(id_val, [1001, 1002, 1003, 1004, 1501, 1502, 1503, 1504])
                fixinj_dose = getdose(em, dose_names{1});
                selected_dose_names = {dose_names{1}};
            else
                fixinj_dose = [getdose(em, dose_names{1}), getdose(em, dose_names{3})];
                selected_dose_names = {dose_names{1}, dose_names{3}};
            end
        else
            fixinj_dose = [getdose(em, dose_names{1}), getdose(em, dose_names{3})];
            selected_dose_names = {dose_names{1}, dose_names{3}};
        end
    elseif strcmp(schedule, 'single dose')
        if any(case_no == [10, 11, 12, 13])
            if ismember(id_val, [1001, 1002, 1003, 1004, 1501, 1502, 1503, 1504])
                fixinj_dose = getdose(em, dose_names{2});
                selected_dose_names = {dose_names{2}};
            else
                fixinj_dose = [getdose(em, dose_names{2}), getdose(em, dose_names{4})];
                selected_dose_names = {dose_names{2}, dose_names{4}};
            end
        else
            fixinj_dose = [getdose(em, dose_names{2}), getdose(em, dose_names{4})];
            selected_dose_names = {dose_names{2}, dose_names{4}};
        end
    elseif strcmp(schedule, 'cont_infusion')
        fixinj_dose = [getdose(em, dose_names{1}), getdose(em, dose_names{2}), getdose(em, dose_names{3})];
        selected_dose_names = {dose_names{1}, dose_names{2}, dose_names{3}};
    else
        error('Unsupported schedule: %s', schedule);
    end

    simData = simulate(em, param_values, fixinj_dose);
    [Tsim, Xsim] = selectbyname(simData, outvec);

    sim_key = sprintf('case%02d_id%d', case_no, id_val);
end
