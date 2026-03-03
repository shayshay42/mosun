function run_phase1_mosu_matlab_reference()
% Phase-1 style mosunetuzumab sweep in MATLAB/SimBiology.
% Uses shared sampled-patient and regimen design from generated/phase1_design.

set(0, 'DefaultFigureVisible', 'off');
warning('off', 'all');

repo_root = fileparts(fileparts(mfilename('fullpath')));
assets_dir = fullfile(repo_root, 'assets');
matlab_dir = fullfile(assets_dir, 'Supp Matlab Code');
design_dir = fullfile(repo_root, 'generated', 'phase1_design');
design_dir_env = getenv('PHASE1_DESIGN_DIR');
if ~isempty(design_dir_env)
    if startsWith(design_dir_env, filesep)
        design_dir = design_dir_env;
    else
        design_dir = fullfile(repo_root, design_dir_env);
    end
end

out_dir = fullfile(repo_root, 'generated', 'phase1_matlab');
out_dir_env = getenv('PHASE1_MATLAB_OUT_DIR');
if ~isempty(out_dir_env)
    if startsWith(out_dir_env, filesep)
        out_dir = out_dir_env;
    else
        out_dir = fullfile(repo_root, out_dir_env);
    end
end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

patients_path = fullfile(design_dir, 'patients.csv');
reg_path = fullfile(design_dir, 'regimen_events.csv');
if ~exist(patients_path, 'file') || ~exist(reg_path, 'file')
    error('Design files missing. Run scripts/generate_phase1_design.py first.');
end

patients_tbl = readtable(patients_path);
reg_tbl = readtable(reg_path);

cd(matlab_dir);
s = sbioloadproject('TDBr26_6_paper.sbproj');
c = struct2cell(s);
model = c{1};

% Reset and activate variants for human DLBCL/tumor setting.
for ii = 1:length(model.variant)
    model.variant(ii).Active = false;
end
variant_ids = [5, 9, 14, 20, 24, 25, 27, 28];
for ii = 1:length(variant_ids)
    model.variant(variant_ids(ii)).Active = true;
end

% Apply DLBCL baseline parameters from the paper parameter table when available.
param_table = readtable(fullfile(assets_dir, 'params_41540_2020_145_MOESM2_ESM.xlsx'), 'Sheet', 'Sheet1');
for ii = 1:height(param_table)
    pname = string(param_table.NAME(ii));
    pval = param_table.DLBCL(ii);
    if ~isnan(pval)
        pobj = sbioselect(model, 'Type', 'parameter', 'Name', char(pname));
        if ~isempty(pobj)
            pobj(1).Value = pval;
        end
    end
end

paramNames = { ...
    'Bpbo_perml', 'Bpbref_perml', 'Trpbo_perml', 'Trpbref_perml', ...
    'KBptumor', 'KTrptumor', 'kBtumorprolif', ...
    'PKflag', 'fvalidation', 'VPid', 'end_time'};
param_objs = SimBiology.Parameter.empty();
for ii = 1:length(paramNames)
    pobj = sbioselect(model, 'Type', 'parameter', 'Name', paramNames{ii});
    if isempty(pobj)
        error('Missing parameter: %s', paramNames{ii});
    end
    param_objs(end+1) = pobj(1); %#ok<AGROW>
end

cs = getconfigset(model, 'active');
cs.StopTime = 84;
cs.MaximumWallClock = 45;

outNames = {'IL6combo', 'Btumor'};
bw_kg_for_dose_conversion = 70.0;

results = cell(height(patients_tbl) * numel(unique(reg_tbl.regimen)), 11);
r = 1;

regimens = unique(reg_tbl.regimen, 'stable');
for rr = 1:numel(regimens)
    reg_name = string(regimens(rr));
    reg_rows = reg_tbl(strcmp(string(reg_tbl.regimen), reg_name), :);
    reg_rows = sortrows(reg_rows, 'event_idx');
    dose_times = reg_rows.time_day';
    dose_mg = reg_rows.dose_mg';
    dose_ugkg = dose_mg * 1000.0 / bw_kg_for_dose_conversion;

    d = sbiodose(sprintf('dose_%s', regexprep(char(reg_name), '[^A-Za-z0-9_]', '_')), 'schedule');
    d.TargetName = 'TDBc_ugperkg';
    d.Time = dose_times;
    d.Amount = dose_ugkg;
    d.Rate = zeros(size(dose_times));

    for pp = 1:height(patients_tbl)
        prow = patients_tbl(pp, :);
        pvals = [ ...
            prow.Bpbo_perml, prow.Bpbref_perml, prow.Trpbo_perml, prow.Trpbref_perml, ...
            prow.KBptumor, prow.KTrptumor, prow.kBtumorprolif, ...
            1, 0, 1, 84];

        try
            for jj = 1:length(param_objs)
                param_objs(jj).Value = pvals(jj);
            end

            simData = sbiosimulate(model, cs, [], d);
            [Tsim, Xsim] = selectbyname(simData, outNames);
            il6 = Xsim(:, 1);
            bt = Xsim(:, 2);
            il6_peak_0_2 = max(il6(Tsim <= 2));
            bt0 = bt(1);
            [Tuniq, Iuniq] = unique(Tsim, 'stable');
            btuniq = bt(Iuniq);
            bt42 = interp1(Tuniq, btuniq, 42.0, 'linear', 'extrap');
            tumor_resid_day42 = bt42 / max(bt0, eps);
            tumor_cfbl_day42 = (bt42 - bt0) / max(bt0, eps);
            loss_raw = 0.5 * log10(1 + max(il6_peak_0_2, 0)) + 0.5 * tumor_resid_day42;
            status = "ok";
        catch ME
            il6_peak_0_2 = NaN;
            bt0 = NaN;
            bt42 = NaN;
            tumor_resid_day42 = NaN;
            tumor_cfbl_day42 = NaN;
            loss_raw = NaN;
            status = "error";
            warning('Simulation failed for regimen=%s patient=%d: %s', reg_name, prow.patient_id, ME.message);
        end

        results(r, :) = { ...
            char(reg_name), char(reg_rows.regimen_type(1)), ...
            double(prow.patient_id), ...
            il6_peak_0_2, bt0, bt42, tumor_resid_day42, tumor_cfbl_day42, ...
            loss_raw, double(prow.BT_ratio_tumor_init), char(status)};
        r = r + 1;
    end
end

results = results(1:r-1, :);
results_tbl = cell2table(results, 'VariableNames', { ...
    'regimen', 'regimen_type', 'patient_id', ...
    'il6_peak_0_2', 'btumor_t0', 'btumor_day42', ...
    'tumor_resid_day42', 'tumor_cfbl_day42', ...
    'loss_raw', 'bt_ratio_tumor_init', 'status'});

out_path = fullfile(out_dir, 'phase1_metrics_matlab.csv');
writetable(results_tbl, out_path);
fprintf('Wrote %s\n', out_path);
end
