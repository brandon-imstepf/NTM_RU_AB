function merge_sweep_outputs(task_dir, out_file)
% MERGE_SWEEP_OUTPUTS  Stitch together per-SLURM-array-task
% sweep_output_task##.mat files (each a NaN-everywhere-except-its-own-slice
% copy of the full sweep grid) into one combined .mat file.
%
% Usage:
%   merge_sweep_outputs()                 % looks in pwd, writes sensitivity_sweep_output.mat
%   merge_sweep_outputs(task_dir)
%   merge_sweep_outputs(task_dir, out_file)

if nargin < 1 || isempty(task_dir), task_dir = pwd; end
if nargin < 2 || isempty(out_file)
    out_file = fullfile(task_dir, 'sensitivity_sweep_output.mat');
end

files = dir(fullfile(task_dir, 'sweep_output_task*.mat'));
assert(~isempty(files), 'No sweep_output_task*.mat files found in %s', task_dir);

first = load(fullfile(files(1).folder, files(1).name));
spread_final = nan(size(first.spread_final));
spread_t50   = nan(size(first.spread_t50));
burden_final = nan(size(first.burden_final));
covered      = false(size(spread_final));

for k = 1:numel(files)
    d = load(fullfile(files(k).folder, files(k).name));
    filled  = ~isnan(d.spread_final);
    overlap = filled & covered;
    assert(~any(overlap(:)), ...
        'Task file %s fills a grid cell another task already filled -- check the SLURM partitioning logic.', ...
        files(k).name);
    spread_final(filled) = d.spread_final(filled);
    spread_t50(filled)   = d.spread_t50(filled);
    burden_final(filled) = d.burden_final(filled);
    covered = covered | filled;
end

fprintf('%d/%d grid cells filled from %d task files.\n', nnz(covered), numel(covered), numel(files));
if ~all(covered(:))
    warning('%d grid cells were never filled by ANY task file -- a task likely crashed or its slice was miscomputed.', ...
        nnz(~covered));
end

coupling_terms = first.coupling_terms; %#ok<NASGU>
kappa_frac     = first.kappa_frac;     %#ok<NASGU>
baseline_val   = first.baseline_val;   %#ok<NASGU>
r_frac         = first.r_frac;         %#ok<NASGU>
Dif_frac       = first.Dif_frac;       %#ok<NASGU>

save(out_file, 'coupling_terms', 'kappa_frac', 'baseline_val', 'r_frac', 'Dif_frac', ...
    'spread_final', 'spread_t50', 'burden_final');
fprintf('Saved merged sweep to %s\n', out_file);

end
