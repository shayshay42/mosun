function export_full_model_reconstruction_bundle()
% Export a complete SimBiology reconstruction bundle for Julia translation.
% Bundle includes:
% - Structural tables: compartments/species/parameters/reactions/rules/events/observables
% - Variants and doses (flattened with all content fields)
% - Sparse stoichiometry with species/reaction names
% - Active configset dump
% - SBML export
% - getequations text for base and selected scenario contexts

set(0, 'DefaultFigureVisible', 'off');
warning('off', 'all');

repo_root = fileparts(fileparts(mfilename('fullpath')));
matlab_dir = fullfile(repo_root, 'assets', 'Supp Matlab Code');
out_dir = fullfile(repo_root, 'generated', 'model_reconstruction_bundle');
eq_dir = fullfile(out_dir, 'equations');

if ~exist(out_dir, 'dir'), mkdir(out_dir); end
if ~exist(eq_dir, 'dir'), mkdir(eq_dir); end

fprintf('Loading model from: %s\n', matlab_dir);
cd(matlab_dir);
addpath(matlab_dir);

s = sbioloadproject('TDBr26_6_paper.sbproj');
c = struct2cell(s);
model = c{1};
cs = getconfigset(model, 'active');

write_summary(model, cs, fullfile(out_dir, 'model_summary.txt'));
write_compartments_table(model, fullfile(out_dir, 'compartments.tsv'));
write_species_table(model, fullfile(out_dir, 'species.tsv'));
write_parameters_table(model, fullfile(out_dir, 'parameters.tsv'));
write_reactions_table(model, fullfile(out_dir, 'reactions.tsv'));
write_rules_table(model, fullfile(out_dir, 'rules.tsv'));
write_events_table(model, fullfile(out_dir, 'events.tsv'));
write_observables_table(model, fullfile(out_dir, 'observables.tsv'));
write_variants_table(model, fullfile(out_dir, 'variants.tsv'));
write_doses_table(model, fullfile(out_dir, 'doses.tsv'));
write_stoich_table(model, fullfile(out_dir, 'stoich_sparse.tsv'));
write_configset_dump(cs, fullfile(out_dir, 'configset_active_dump.txt'));

sbml_path = fullfile(out_dir, 'TDBr26_6_paper.sbml');
try
    sbmlexport(model, sbml_path);
    fprintf('Wrote %s\n', sbml_path);
catch ME
    warning('SBML export failed: %s', ME.message);
end

write_equations_bundle(model, cs, eq_dir, fullfile(out_dir, 'equation_manifest.tsv'));

write_text(fullfile(out_dir, 'README.txt'), bundle_readme_text());

fprintf('\nDone. Reconstruction bundle written to:\n%s\n', out_dir);
end


function write_summary(model, cs, out_path)
f = fopen(out_path, 'w');
fprintf(f, 'MODEL_NAME\t%s\n', as_char(model.Name));
fprintf(f, 'COUNTS\tcompartments\t%d\tspecies\t%d\tparameters\t%d\treactions\t%d\trules\t%d\tevents\t%d\tobservables\t%d\tvariants\t%d\tdoses\t%d\n', ...
    length(model.Compartments), length(model.Species), length(model.Parameters), ...
    length(model.Reactions), length(model.Rules), length(model.Events), ...
    length(model.Observables), length(model.Variants), length(model.Doses));
fprintf(f, 'CONFIGSET\tname\t%s\tsolver\t%s\tstop_time\t%g\tmax_wall_clock\t%g\n', ...
    as_char(cs.Name), as_char(cs.SolverType), double(cs.StopTime), double(cs.MaximumWallClock));
fclose(f);
fprintf('Wrote %s\n', out_path);
end


function write_compartments_table(model, out_path)
n = length(model.Compartments);
rows = cell(n, 9);
for i = 1:n
    cp = model.Compartments(i);
    rows{i,1} = i;
    rows{i,2} = as_char(cp.Name);
    rows{i,3} = as_char(cp.Parent.Name);
    rows{i,4} = double(cp.Capacity);
    rows{i,5} = as_char(cp.CapacityUnits);
    rows{i,6} = double(cp.ConstantCapacity);
    rows{i,7} = as_char(cp.Type);
    rows{i,8} = as_char(cp.Notes);
    rows{i,9} = as_char(cp.Tag);
end
tbl = cell2table(rows, 'VariableNames', { ...
    'idx', 'name', 'parent', 'capacity', 'capacity_units', ...
    'constant_capacity', 'type', 'notes', 'tag'});
writetable_tsv(tbl, out_path);
end


function write_species_table(model, out_path)
n = length(model.Species);
rows = cell(n, 12);
for i = 1:n
    sp = model.Species(i);
    rows{i,1} = i;
    rows{i,2} = as_char(sp.Name);
    rows{i,3} = as_char(sp.Parent.Name);
    rows{i,4} = double(sp.InitialAmount);
    rows{i,5} = as_char(sp.InitialAmountUnits);
    rows{i,6} = as_char(sp.Units);
    rows{i,7} = double(sp.ConstantAmount);
    rows{i,8} = double(sp.BoundaryCondition);
    rows{i,9} = double(sp.Constant);
    rows{i,10} = as_char(sp.Type);
    rows{i,11} = as_char(sp.Notes);
    rows{i,12} = as_char(sp.Tag);
end
tbl = cell2table(rows, 'VariableNames', { ...
    'idx', 'name', 'parent', 'initial_amount', 'initial_amount_units', ...
    'units', 'constant_amount', 'boundary_condition', 'constant', ...
    'type', 'notes', 'tag'});
writetable_tsv(tbl, out_path);
end


function write_parameters_table(model, out_path)
n = length(model.Parameters);
rows = cell(n, 11);
for i = 1:n
    p = model.Parameters(i);
    rows{i,1} = i;
    rows{i,2} = as_char(p.Name);
    rows{i,3} = as_char(p.Parent.Name);
    rows{i,4} = double(p.Value);
    rows{i,5} = as_char(p.ValueUnits);
    rows{i,6} = as_char(p.Units);
    rows{i,7} = double(p.ConstantValue);
    rows{i,8} = double(p.BoundaryCondition);
    rows{i,9} = double(p.Constant);
    rows{i,10} = as_char(p.Notes);
    rows{i,11} = as_char(p.Tag);
end
tbl = cell2table(rows, 'VariableNames', { ...
    'idx', 'name', 'parent', 'value', 'value_units', 'units', ...
    'constant_value', 'boundary_condition', 'constant', 'notes', 'tag'});
writetable_tsv(tbl, out_path);
end


function write_reactions_table(model, out_path)
n = length(model.Reactions);
rows = cell(n, 10);
for i = 1:n
    rxn = model.Reactions(i);
    rows{i,1} = i;
    rows{i,2} = as_char(rxn.Name);
    rows{i,3} = as_char(rxn.Reaction);
    rows{i,4} = as_char(rxn.ReactionRate);
    rows{i,5} = double(rxn.Reversible);
    rows{i,6} = double(rxn.Active);
    rows{i,7} = class_or_empty(rxn.KineticLaw);
    rows{i,8} = json_or_empty(param_var_names_or_empty(rxn.KineticLaw));
    rows{i,9} = json_or_empty(kinetic_law_params_or_empty(rxn.KineticLaw));
    rows{i,10} = as_char(rxn.Tag);
end
tbl = cell2table(rows, 'VariableNames', { ...
    'idx', 'name', 'reaction', 'reaction_rate', 'reversible', 'active', ...
    'kinetic_law_class', 'kinetic_law_param_names_json', ...
    'kinetic_law_param_values_json', 'tag'});
writetable_tsv(tbl, out_path);
end


function write_rules_table(model, out_path)
n = length(model.Rules);
rows = cell(n, 6);
for i = 1:n
    rr = model.Rules(i);
    rows{i,1} = i;
    rows{i,2} = as_char(rr.Name);
    rows{i,3} = as_char(rr.RuleType);
    rows{i,4} = as_char(rr.Rule);
    rows{i,5} = double(rr.Active);
    rows{i,6} = as_char(rr.Tag);
end
tbl = cell2table(rows, 'VariableNames', {'idx', 'name', 'rule_type', 'rule', 'active', 'tag'});
writetable_tsv(tbl, out_path);
end


function write_events_table(model, out_path)
n = length(model.Events);
rows = cell(n, 7);
for i = 1:n
    ev = model.Events(i);
    rows{i,1} = i;
    rows{i,2} = as_char(ev.Name);
    rows{i,3} = as_char(ev.Trigger);
    rows{i,4} = json_or_empty(ev.EventFcns);
    rows{i,5} = double(ev.Active);
    rows{i,6} = as_char(ev.Type);
    rows{i,7} = as_char(ev.Tag);
end
tbl = cell2table(rows, 'VariableNames', {'idx', 'name', 'trigger', 'event_fcns_json', 'active', 'type', 'tag'});
writetable_tsv(tbl, out_path);
end


function write_observables_table(model, out_path)
n = length(model.Observables);
rows = cell(n, 6);
for i = 1:n
    ob = model.Observables(i);
    rows{i,1} = i;
    rows{i,2} = as_char(ob.Name);
    rows{i,3} = as_char(ob.Expression);
    rows{i,4} = double(ob.Active);
    rows{i,5} = as_char(ob.Type);
    rows{i,6} = as_char(ob.Tag);
end
tbl = cell2table(rows, 'VariableNames', {'idx', 'name', 'expression', 'active', 'type', 'tag'});
writetable_tsv(tbl, out_path);
end


function write_variants_table(model, out_path)
rows = {};
r = 1;
for i = 1:length(model.Variants)
    v = model.Variants(i);
    C = v.Content;
    if isempty(C)
        rows(r,:) = {i, as_char(v.Name), 0, '', '', '', NaN, '', 0}; %#ok<AGROW>
        r = r + 1;
        continue;
    end

    for j = 1:size(C,1)
        entry = C{j};
        action = as_char(entry{1});
        entry_class = as_char(entry{2});
        entry_name = as_char(entry{3});
        entry_value = entry{4};
        [vnum, is_num] = scalar_to_double(entry_value);
        rows(r,:) = { ...
            i, as_char(v.Name), j, action, entry_class, entry_name, ...
            vnum, value_to_text(entry_value), is_num}; %#ok<AGROW>
        r = r + 1;
    end
end
tbl = cell2table(rows, 'VariableNames', { ...
    'variant_idx', 'variant_name', 'entry_idx', 'action', 'class', ...
    'name', 'value_numeric', 'value_text', 'is_numeric'});
writetable_tsv(tbl, out_path);
end


function write_doses_table(model, out_path)
n = length(model.Doses);
rows = cell(n, 16);
for i = 1:n
    d = model.Doses(i);
    rows{i,1} = i;
    rows{i,2} = as_char(d.Name);
    rows{i,3} = class(d);
    rows{i,4} = as_char(d.TargetName);
    rows{i,5} = json_or_empty(get_prop_or_default(d, 'Amount', []));
    rows{i,6} = json_or_empty(get_prop_or_default(d, 'Rate', []));
    rows{i,7} = json_or_empty(get_prop_or_default(d, 'Time', []));
    rows{i,8} = json_or_empty(get_prop_or_default(d, 'StartTime', []));
    rows{i,9} = json_or_empty(get_prop_or_default(d, 'Interval', []));
    rows{i,10} = json_or_empty(get_prop_or_default(d, 'RepeatCount', []));
    rows{i,11} = as_char(get_prop_or_default(d, 'AmountUnits', ''));
    rows{i,12} = as_char(get_prop_or_default(d, 'RateUnits', ''));
    rows{i,13} = as_char(get_prop_or_default(d, 'TimeUnits', ''));
    rows{i,14} = as_char(get_prop_or_default(d, 'EventMode', ''));
    rows{i,15} = double(get_prop_or_default(d, 'Active', 1));
    rows{i,16} = as_char(get_prop_or_default(d, 'Tag', ''));
end
tbl = cell2table(rows, 'VariableNames', { ...
    'idx', 'name', 'class', 'target_name', 'amount_json', 'rate_json', 'time_json', ...
    'start_time_json', 'interval_json', 'repeat_count_json', ...
    'amount_units', 'rate_units', 'time_units', 'event_mode', 'active', 'tag'});
writetable_tsv(tbl, out_path);
end


function write_stoich_table(model, out_path)
S = getstoichmatrix(model);
[ri, ci, vi] = find(S);
n = numel(ri);
rows = cell(n, 5);
for i = 1:n
    rows{i,1} = ri(i);
    rows{i,2} = ci(i);
    rows{i,3} = vi(i);
    rows{i,4} = as_char(model.Species(ri(i)).Name);
    rows{i,5} = as_char(model.Reactions(ci(i)).Name);
end
tbl = cell2table(rows, 'VariableNames', { ...
    'species_idx', 'reaction_idx', 'stoich_coeff', 'species_name', 'reaction_name'});
writetable_tsv(tbl, out_path);
end


function write_configset_dump(cs, out_path)
txt = evalc('disp(cs)');
write_text(out_path, txt);
fprintf('Wrote %s\n', out_path);
end


function write_equations_bundle(model, cs, eq_dir, manifest_path)
rows = {};
r = 1;

    function add_scenario(scenario_name, variant_ids, doses_obj, notes)
        safe_name = sanitize_filename(scenario_name);
        eq_path = fullfile(eq_dir, [safe_name, '.txt']);
        if isempty(variant_ids)
            vv = [];
            vnames = {};
        else
            vv = model.Variants(variant_ids);
            vnames = arrayfun(@(x) as_char(x.Name), vv, 'UniformOutput', false);
        end

        if isempty(doses_obj)
            dd = [];
            dnames = {};
        else
            dd = doses_obj;
            dnames = arrayfun(@(x) as_char(x.Name), dd, 'UniformOutput', false);
        end

        eq_txt = string(getequations(model, cs, vv, dd));
        write_text(eq_path, char(eq_txt));

        rows(r,:) = { ...
            scenario_name, json_or_empty(variant_ids), json_or_empty(vnames), ...
            json_or_empty(dnames), notes, eq_path}; %#ok<AGROW>
        r = r + 1;
    end

% Base equations (no variants/doses applied).
add_scenario('base_no_variants_no_doses', [], [], ...
    'Raw model equations using active configset only.');

% Legacy paper contexts.
add_scenario('legacy_case10_variants_dosepool', [5 9 11 13 14 24], model.Doses([7 8 18 19]), ...
    'Variants+dose pool used in MainRun case 10-13.');
add_scenario('legacy_case20_variants_dosepool', [5 9 11 13 14 15 24], model.Doses([7 8 18 19]), ...
    'Variants+dose pool used in MainRun case 20-22.');
add_scenario('legacy_case30_variants_dosepool', [5 9 11 13 14 20 22 24], model.Doses([8 19 20]), ...
    'Variants+dose pool used in MainRun case 30.');

% Human phase-1 variant context without/with representative schedules.
phase1_vid = [5 9 14 20 24 25 27 28];
add_scenario('phase1_human_variants_no_dose', phase1_vid, [], ...
    'Variants used by phase1 scripts, no exogenous dose.');

d_fixed = sbiodose('phase1_fixed_0p8mg_q3w', 'schedule');
d_fixed.TargetName = 'TDBc_ugperkg';
d_fixed.Time = [0 21 42 63];
d_fixed.Amount = (0.8 * 1000 / 70) * ones(size(d_fixed.Time));
d_fixed.Rate = zeros(size(d_fixed.Time));
add_scenario('phase1_human_fixed_0p8mg_q3w', phase1_vid, d_fixed, ...
    'Representative fixed-q3w schedule for equation expansion.');

d_step = sbiodose('phase1_step_0p8_2_6mg', 'schedule');
d_step.TargetName = 'TDBc_ugperkg';
d_step.Time = [0 7 14 21 42 63];
d_step.Amount = ([0.8 2 6 6 6 6] * 1000 / 70);
d_step.Rate = zeros(size(d_step.Time));
add_scenario('phase1_human_step_0p8_2_6mg_q3w', phase1_vid, d_step, ...
    'Representative cycle-1 step-up + q3w target schedule.');

% One equation export per model-defined dose in the project.
for i = 1:length(model.Doses)
    dose_name = as_char(model.Doses(i).Name);
    sc_name = sprintf('base_modeldose_%02d_%s', i, dose_name);
    note = sprintf('Base model with only project dose idx=%d active.', i);
    add_scenario(sc_name, [], model.Doses(i), note);
end

tbl = cell2table(rows, 'VariableNames', { ...
    'scenario_name', 'variant_ids_json', 'variant_names_json', ...
    'dose_names_json', 'notes', 'equation_file'});
writetable_tsv(tbl, manifest_path);
fprintf('Wrote %s\n', manifest_path);
end


function txt = bundle_readme_text()
txt = [ ...
    "This folder is a complete SimBiology export bundle for Julia reconstruction.", newline, ...
    newline, ...
    "Key files:", newline, ...
    "- model_summary.txt: counts and active configset summary", newline, ...
    "- compartments/species/parameters/reactions/rules/events/observables TSVs", newline, ...
    "- variants.tsv and doses.tsv: flattened variant content and dose metadata", newline, ...
    "- stoich_sparse.tsv: sparse stoichiometry entries with indices and names", newline, ...
    "- configset_active_dump.txt: active SimBiology configset dump", newline, ...
    "- TDBr26_6_paper.sbml: SBML export (if sbmlexport succeeds)", newline, ...
    "- equations/*.txt + equation_manifest.tsv: getequations text per scenario", newline, ...
    newline, ...
    "All dose amounts in phase1 representative scenarios are converted mg->ug/kg using BW=70 kg.", newline];
end


function names = param_var_names_or_empty(kinetic_law)
names = {};
if isempty(kinetic_law)
    return;
end
try
    names = kinetic_law.ParameterVariableNames;
catch
    names = {};
end
end


function entries = kinetic_law_params_or_empty(kinetic_law)
entries = struct('name', {}, 'value', {});
if isempty(kinetic_law)
    return;
end
try
    kp = kinetic_law.Parameters;
    entries = repmat(struct('name', '', 'value', NaN), length(kp), 1);
    for i = 1:length(kp)
        entries(i).name = as_char(kp(i).Name);
        entries(i).value = double(kp(i).Value);
    end
catch
    entries = struct('name', {}, 'value', {});
end
end


function out = class_or_empty(x)
if isempty(x)
    out = '';
else
    out = class(x);
end
end


function out = get_prop_or_default(obj, prop_name, default_val)
if isprop(obj, prop_name)
    out = obj.(prop_name);
else
    out = default_val;
end
end


function [vnum, is_num] = scalar_to_double(v)
is_num = 0;
vnum = NaN;
if isnumeric(v) && isscalar(v)
    is_num = 1;
    vnum = double(v);
    return;
end
if islogical(v) && isscalar(v)
    is_num = 1;
    vnum = double(v);
    return;
end
if ischar(v) || isstring(v)
    tmp = str2double(char(v));
    if ~isnan(tmp)
        is_num = 1;
        vnum = tmp;
    end
end
end


function out = value_to_text(v)
if ischar(v)
    out = v;
elseif isstring(v)
    out = char(v);
elseif isnumeric(v) || islogical(v)
    out = jsonencode(v);
else
    out = char(string(v));
end
end


function out = json_or_empty(v)
try
    if isempty(v)
        out = '[]';
    else
        out = char(jsonencode(v));
    end
catch
    out = '[]';
end
end


function out = as_char(x)
if isstring(x)
    out = char(x);
elseif ischar(x)
    out = x;
elseif isempty(x)
    out = '';
else
    out = char(string(x));
end
end


function s = sanitize_filename(name)
s = as_char(name);
s = regexprep(s, '[^A-Za-z0-9_]+', '_');
s = regexprep(s, '_+', '_');
s = regexprep(s, '^_+', '');
s = regexprep(s, '_+$', '');
if isempty(s)
    s = 'unnamed';
end
end


function write_text(path, txt)
f = fopen(path, 'w');
if f < 0
    error('Could not open for writing: %s', path);
end
fprintf(f, '%s', txt);
fclose(f);
fprintf('Wrote %s\n', path);
end


function writetable_tsv(tbl, out_path)
writetable(tbl, out_path, 'FileType', 'text', 'Delimiter', '\t');
fprintf('Wrote %s\n', out_path);
end
