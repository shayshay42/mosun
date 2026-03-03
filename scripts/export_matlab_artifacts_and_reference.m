set(0, 'DefaultFigureVisible', 'off');

repo_root = fileparts(fileparts(mfilename('fullpath')));
matlab_dir = fullfile(repo_root, 'assets', 'Supp Matlab Code');
out_root = fullfile(repo_root, 'generated');
out_ref_dir = fullfile(out_root, 'matlab_reference');
model_json_path = fullfile(out_root, 'model_export.json');
manifest_path = fullfile(out_ref_dir, 'manifest.csv');

if ~exist(out_root, 'dir'), mkdir(out_root); end
if ~exist(out_ref_dir, 'dir'), mkdir(out_ref_dir); end

fprintf('Working directory: %s\n', matlab_dir);
cd(matlab_dir);

s = sbioloadproject('TDBr26_6_paper.sbproj');
c = struct2cell(s);
model = c{1};

export_model_json(model, model_json_path);
fprintf('Wrote model artifact: %s\n', model_json_path);

cases = [10, 11, 12, 13, 20, 21, 22, 30];
meta_rows = {};
row = 1;

for case_no = cases
    fprintf('\n=== CASE %d ===\n', case_no);

    [variant_id, dose_id, outvec, filename, paramNames, idvec] = case_config(case_no);

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

    for idpos = 1:length(idvec)
        this_id = idvec(idpos);
        pos = find(data.ID == this_id);
        if isempty(pos)
            warning('Skipping case %d ID %d: no rows in %s', case_no, this_id, filename);
            continue;
        end

        Tdata_tmp = sort(unique(data.Time(pos)));
        SimTime = max(Tdata_tmp) + 0.1;
        if case_no == 30
            SimTime = max(Tdata_tmp);
        end

        outtime = sort(unique([(0:0.1:SimTime)'; Tdata_tmp]));
        outtime = outtime(outtime >= 0);

        em = export(model, parameters);
        em.SimulationOptions.OutputTimes = outtime;
        em.SimulationOptions.StopTime = SimTime;
        em.SimulationOptions.MaximumWallClock = 120;

        init_Bpbo = mean(data.Bcells(pos(data.Time(pos)==0)), 'omitnan') * 1000;
        init_Tpbo = mean(data.CD8Tcells(pos(data.Time(pos)==0)), 'omitnan') * 1000;

        if case_no ~= 30
            param_values = [init_Bpbo, init_Bpbo, init_Tpbo, init_Tpbo, 0, this_id, SimTime];
        else
            param_values = [init_Bpbo, init_Bpbo, init_Tpbo, init_Tpbo, 1];
        end

        schedule = schedule_for_id(data.schedule(pos));
        [fixinj_dose, selected_dose_names] = select_dose_for_case(case_no, this_id, schedule, em, dose_names);

        simData = simulate(em, param_values, fixinj_dose);
        [Tsim, Xsim] = selectbyname(simData, outvec);

        sim_key = sprintf('case%02d_id%d', case_no, this_id);
        out_csv = fullfile(out_ref_dir, [sim_key, '.csv']);
        out_tbl = array2table([Tsim, Xsim], 'VariableNames', ['time', outvec']);
        writetable(out_tbl, out_csv);

        rel_csv = fullfile('generated', 'matlab_reference', [sim_key, '.csv']);
        meta_rows{row,1} = sim_key;
        meta_rows{row,2} = case_no;
        meta_rows{row,3} = this_id;
        meta_rows{row,4} = schedule;
        meta_rows{row,5} = SimTime;
        meta_rows{row,6} = jsonencode(outvec);
        meta_rows{row,7} = jsonencode(variant_id);
        meta_rows{row,8} = jsonencode(dose_names);
        meta_rows{row,9} = jsonencode(selected_dose_names);
        meta_rows{row,10} = jsonencode(paramNames);
        meta_rows{row,11} = jsonencode(param_values);
        meta_rows{row,12} = rel_csv;
        row = row + 1;

        fprintf('Saved %-16s | ntime=%4d | schedule=%s\n', sim_key, length(Tsim), schedule);
    end

    for ii = 1:length(variant_id)
        model.variant(variant_id(ii)).Active = false;
    end
    clear parameters;
end

meta_tbl = cell2table(meta_rows, 'VariableNames', {
    'sim_key', 'case_no', 'id', 'schedule', 'sim_time', 'outvec_json', ...
    'variant_ids_json', 'dose_pool_json', 'selected_doses_json', ...
    'param_names_json', 'param_values_json', 'sim_csv_relpath'});
writetable(meta_tbl, manifest_path);

fprintf('\nWrote manifest: %s\n', manifest_path);
fprintf('Total simulations exported: %d\n', height(meta_tbl));


function [variant_id, dose_id, outvec, filename, paramNames, idvec] = case_config(case_no)
    switch case_no
        case {10, 11, 12, 13}
            variant_id = [5; 9; 11; 13; 14; 24];
            dose_id = [7, 8, 18, 19];
            outvec = {'Bpb_perml', 'totTpb_perml', 'Tafraction_pb', 'BTtotRatio_tiss', 'BTtotRatio_tiss2', 'IL6combo'};
            filename = 'dataset_1.xlsx';
            paramNames = {'Bpbo_perml', 'Bpbref_perml', 'Trpbo_perml', 'Trpbref_perml', 'PKflag', 'VPid', 'end_time'};
        case {20, 21, 22}
            variant_id = [5; 9; 11; 13; 14; 15; 24];
            dose_id = [7, 8, 18, 19];
            outvec = {'Bpb_perml', 'totTpb_perml', 'Tafraction_pb'};
            filename = 'dataset_2.xlsx';
            paramNames = {'Bpbo_perml', 'Bpbref_perml', 'Trpbo_perml', 'Trpbref_perml', 'PKflag', 'VPid', 'end_time'};
        case {30}
            variant_id = [5; 9; 11; 13; 14; 20; 22; 24];
            dose_id = [8, 19, 20];
            outvec = {'Bpb_perml', 'totTpb_perml', 'Tafraction_pb'};
            filename = 'dataset_blin.xlsx';
            paramNames = {'Bpbo_perml', 'Bpbref_perml', 'Trpbo_perml', 'Trpbref_perml', 'PKflag'};
        otherwise
            error('Unsupported case number: %d', case_no);
    end

    switch case_no
        case 10
            idvec = [1001, 1002, 1003, 1004, 1501, 1502, 1503, 1504];
        case 11
            idvec = [2001, 2002, 2501, 2502];
        case 12
            idvec = [3001, 3002, 3501, 3502];
        case 13
            idvec = [4001, 4002, 4003, 4004, 4501, 4502, 4503, 20013, 20023, 20033, 20043];
        case 20
            idvec = [3001, 3002, 3003];
        case 21
            idvec = [2001, 2002, 2003, 2004];
        case 22
            idvec = [1001, 1002, 1003, 1004];
        case 30
            idvec = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20];
        otherwise
            idvec = [];
    end
end


function schedule = schedule_for_id(schedule_values)
    u = unique(schedule_values);
    if iscell(u)
        schedule = u{1};
    else
        schedule = u;
    end
    schedule = char(schedule);
end


function [fixinj_dose, selected_dose_names] = select_dose_for_case(case_no, id_val, schedule, em, dose_names)
    switch schedule
        case '4qw'
            switch case_no
                case {10, 11, 12, 13}
                    if ismember(id_val, [1001, 1002, 1003, 1004, 1501, 1502, 1503, 1504])
                        fixinj_dose = getdose(em, dose_names{1});
                        selected_dose_names = {dose_names{1}};
                    else
                        fixinj_dose = [getdose(em, dose_names{1}), getdose(em, dose_names{3})];
                        selected_dose_names = {dose_names{1}, dose_names{3}};
                    end
                case {20, 21, 22}
                    fixinj_dose = [getdose(em, dose_names{1}), getdose(em, dose_names{3})];
                    selected_dose_names = {dose_names{1}, dose_names{3}};
                otherwise
                    fixinj_dose = getdose(em, dose_names{1});
                    selected_dose_names = {dose_names{1}};
            end
        case 'single dose'
            switch case_no
                case {10, 11, 12, 13}
                    if ismember(id_val, [1001, 1002, 1003, 1004, 1501, 1502, 1503, 1504])
                        fixinj_dose = getdose(em, dose_names{2});
                        selected_dose_names = {dose_names{2}};
                    else
                        fixinj_dose = [getdose(em, dose_names{2}), getdose(em, dose_names{4})];
                        selected_dose_names = {dose_names{2}, dose_names{4}};
                    end
                case {20, 21, 22}
                    fixinj_dose = [getdose(em, dose_names{2}), getdose(em, dose_names{4})];
                    selected_dose_names = {dose_names{2}, dose_names{4}};
                otherwise
                    fixinj_dose = getdose(em, dose_names{2});
                    selected_dose_names = {dose_names{2}};
            end
        case 'cont_infusion'
            fixinj_dose = [getdose(em, dose_names{1}), getdose(em, dose_names{2}), getdose(em, dose_names{3})];
            selected_dose_names = {dose_names{1}, dose_names{2}, dose_names{3}};
        otherwise
            error('Unexpected schedule string: %s', schedule);
    end
end


function export_model_json(model, model_json_path)
    artifact = struct();

    species = repmat(struct('name', '', 'initial_amount', 0), length(model.Species), 1);
    for i = 1:length(model.Species)
        species(i).name = model.Species(i).Name;
        species(i).initial_amount = model.Species(i).InitialAmount;
    end
    artifact.species = species;

    parameters = repmat(struct('name', '', 'value', 0), length(model.Parameters), 1);
    for i = 1:length(model.Parameters)
        parameters(i).name = model.Parameters(i).Name;
        parameters(i).value = model.Parameters(i).Value;
    end
    artifact.parameters = parameters;

    rules = repmat(struct('type', '', 'rule', ''), length(model.Rules), 1);
    for i = 1:length(model.Rules)
        rules(i).type = model.Rules(i).RuleType;
        rules(i).rule = model.Rules(i).Rule;
    end
    artifact.rules = rules;

    reactions = repmat(struct('reaction', '', 'rate', ''), length(model.Reactions), 1);
    for i = 1:length(model.Reactions)
        reactions(i).reaction = model.Reactions(i).Reaction;
        reactions(i).rate = model.Reactions(i).ReactionRate;
    end
    artifact.reactions = reactions;

    S = getstoichmatrix(model);
    [r, c, v] = find(S);
    stoich = struct();
    stoich.rows = r;
    stoich.cols = c;
    stoich.vals = v;
    stoich.n_species = size(S, 1);
    stoich.n_reactions = size(S, 2);
    artifact.stoich = stoich;

    variants = repmat(struct('name', '', 'entries', []), length(model.Variants), 1);
    for i = 1:length(model.Variants)
        variants(i).name = model.Variants(i).Name;
        C = model.Variants(i).Content;
        entries = repmat(struct('action', '', 'class', '', 'name', '', 'value', []), size(C,1), 1);
        for j = 1:size(C,1)
            entry = C{j};
            entries(j).action = as_char(entry{1});
            entries(j).class = as_char(entry{2});
            entries(j).name = as_char(entry{3});
            entries(j).value = entry{4};
        end
        variants(i).entries = entries;
    end
    artifact.variants = variants;

    doses = repmat(struct('name', '', 'class', '', 'target', '', 'amount', [], 'rate', [], ...
        'interval', [], 'repeatcount', [], 'starttime', [], 'time', []), length(model.Doses), 1);
    for i = 1:length(model.Doses)
        d = model.Doses(i);
        doses(i).name = d.Name;
        doses(i).class = class(d);
        doses(i).target = as_char(d.TargetName);
        doses(i).amount = d.Amount;
        doses(i).rate = d.Rate;
        if isprop(d, 'Interval')
            doses(i).interval = d.Interval;
            doses(i).repeatcount = d.RepeatCount;
            doses(i).starttime = d.StartTime;
        end
        if isprop(d, 'Time')
            doses(i).time = d.Time;
        end
    end
    artifact.doses = doses;

    fid = fopen(model_json_path, 'w');
    fprintf(fid, '%s', jsonencode(artifact));
    fclose(fid);
end


function out = as_char(x)
    if isstring(x)
        out = char(x);
    elseif ischar(x)
        out = x;
    else
        out = char(string(x));
    end
end
