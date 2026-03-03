function export_phase1_standard_traces_matlab()
% Export MATLAB trajectories for phase-1 regimens (toxicity/tumor proxies).
% Uses PHASE1_DESIGN_DIR and PHASE1_MATLAB_TRACE_OUT_DIR when provided.

set(0, 'DefaultFigureVisible', 'off');
warning('off', 'all');

repo_root = fileparts(fileparts(mfilename('fullpath')));
assets_dir = fullfile(repo_root, 'assets');
matlab_dir = fullfile(assets_dir, 'Supp Matlab Code');

design_dir = fullfile(repo_root, 'generated', 'phase1_design_standard');
design_dir_env = getenv('PHASE1_DESIGN_DIR');
if ~isempty(design_dir_env)
    if startsWith(design_dir_env, filesep)
        design_dir = design_dir_env;
    else
        design_dir = fullfile(repo_root, design_dir_env);
    end
end

out_dir = fullfile(repo_root, 'generated', 'phase1_standard_matlab_traces');
out_dir_env = getenv('PHASE1_MATLAB_TRACE_OUT_DIR');
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
    error('Design files missing in %s', design_dir);
end

patients_tbl = readtable(patients_path);
reg_tbl = readtable(reg_path);

cd(matlab_dir);
s = sbioloadproject('TDBr26_6_paper.sbproj');
c = struct2cell(s);
model = c{1};

for ii = 1:length(model.variant)
    model.variant(ii).Active = false;
end
variant_ids = [5, 9, 14, 20, 24, 25, 27, 28];
for ii = 1:length(variant_ids)
    model.variant(variant_ids(ii)).Active = true;
end

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
tgrid = unique([0:0.01:2, 2:0.1:84])';

manifest = cell(height(patients_tbl) * numel(unique(reg_tbl.regimen)), 4);
k = 1;

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

        safe_reg = regexprep(char(reg_name), '[^A-Za-z0-9_]', '_');
        fname = sprintf('regimen_%s_patient_%d.csv', safe_reg, double(prow.patient_id));
        out_path = fullfile(out_dir, fname);

        try
            for jj = 1:length(param_objs)
                param_objs(jj).Value = pvals(jj);
            end

            simData = sbiosimulate(model, cs, [], d);
            [Tsim, Xsim] = selectbyname(simData, outNames);
            [Tuniq, Iuniq] = unique(Tsim, 'stable');
            il6u = Xsim(Iuniq, 1);
            btu = Xsim(Iuniq, 2);
            il6g = interp1(Tuniq, il6u, tgrid, 'linear', 'extrap');
            btg = interp1(Tuniq, btu, tgrid, 'linear', 'extrap');

            out_tbl = table(tgrid, il6g, btg, ...
                repmat(string(reg_name), numel(tgrid), 1), ...
                repmat(double(prow.patient_id), numel(tgrid), 1), ...
                'VariableNames', {'time', 'IL6combo', 'Btumor', 'regimen', 'patient_id'});
            writetable(out_tbl, out_path);
            status = "ok";
        catch
            status = "error";
        end

        manifest(k, :) = {char(reg_name), double(prow.patient_id), out_path, char(status)};
        k = k + 1;
    end
end

manifest = manifest(1:k-1, :);
manifest_tbl = cell2table(manifest, 'VariableNames', {'regimen', 'patient_id', 'trace_csv', 'status'});
writetable(manifest_tbl, fullfile(out_dir, 'manifest.csv'));
fprintf('Wrote %s\n', fullfile(out_dir, 'manifest.csv'));
end

