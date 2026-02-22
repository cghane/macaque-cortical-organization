%% performs NaN-safe PCA on voxelwise receptor maps from RM_dump_nan.1D, saves voxelwise PC score maps as 1D files, and stores PCA loadings and metadata.
%% ================================================================
%   neuro_PCA_1D_maps_nanSafe.m
%
%   - Load RM_dump_nan.1D (XYZ + 14 receptor columns)
%   - Robust PCA on receptor columns (4–17) with NaN handling
%   - Save PC scores as 1D files (voxel-wise PC maps)
%   - Output directory:
%       /mnt/scratch/NHP4CYRUS/output/neuro_PCA_maps/
% ================================================================

clear; clc;

%% =========================
%  Paths & I/O
% ==========================
prefix   = '/mnt/scratch/NHP4CYRUS/';
rm_file  = fullfile(prefix, 'RM_dump_nan.1D');

% Create a sensible output directory under "output"
output_dir = ensure_outdir(fullfile(prefix, 'output', 'neuro_PCA_maps'));
fprintf('PC maps will be saved in: %s\n', output_dir);

%% =========================
%  Load neurotransmitter data
% ==========================
if ~isfile(rm_file)
    error('File not found: %s', rm_file);
end

neuro_data = load(rm_file);   % Expecting 17 columns: XYZ(1:3) + receptors(4:17)

if size(neuro_data, 2) ~= 17
    error('RM_dump_nan.1D should have 17 columns (got %d).', size(neuro_data, 2));
end

xyz_vox = neuro_data(:,1:3);   %#ok<NASGU> % keep for potential later use
X_raw   = neuro_data(:,4:17);  % 14 receptor columns

% Receptor names corresponding to cols 4–17
receptor_names_all = {'AMPA_Glutamate', 'Kainate_Glutamate', 'NMDA_Glutamate', ...
                      'GABA_A', 'GABA_A_BZ', 'GABA_B', ...
                      'M1_Acetylcholine', 'M2_Acetylcholine', 'M3_Acetylcholine', ...
                      'HT1A_Serotonin', 'HT2A_Serotonin', ...
                      'Alpha1_Noradrenaline', 'Alpha2_Noradrenaline', ...
                      'D1_Dopamine'};

%% =========================
%  Inspect NaNs & drop bad receptors if needed
% ==========================
nVox = size(X_raw,1);
nan_counts = sum(isnan(X_raw), 1);
finite_counts = sum(isfinite(X_raw), 1);

fprintf('--- Receptor data overview ---\n');
for j = 1:numel(receptor_names_all)
    fprintf('%2d: %-20s  finite=%6d  NaN=%6d\n', ...
        j, receptor_names_all{j}, finite_counts(j), nan_counts(j));
end

% Drop receptors with almost no valid data (e.g. < 10 finite voxels)
min_valid_vox = 10;
valid_cols = finite_counts >= min_valid_vox;

if ~any(valid_cols)
    error('No receptor has at least %d valid voxels. Cannot run PCA.', min_valid_vox);
end

if any(~valid_cols)
    fprintf('\nDropping receptors with < %d valid voxels:\n', min_valid_vox);
    disp(receptor_names_all(~valid_cols)');
end

X_use          = X_raw(:, valid_cols);
receptor_names = receptor_names_all(valid_cols);
nRec           = numel(receptor_names);

fprintf('\nUsing %d / %d receptors for PCA.\n', nRec, numel(receptor_names_all));

%% =========================
%  Z-score with NaN-safe handling
% ==========================
% Compute mean/std ignoring NaNs
mu    = mean(X_use, 1, 'omitnan');
sigma = std(X_use, 0, 1, 'omitnan');

% Avoid division by zero (constant columns)
sigma(sigma == 0 | isnan(sigma)) = 1;

% Center & scale
Xz = (X_use - mu) ./ sigma;

% At this point Xz still has NaNs where original data were NaN.
% Replace NaNs with 0 = "at-mean" after z-scoring.
Xz(isnan(Xz)) = 0;

% Also guard against infs just in case
Xz(~isfinite(Xz)) = 0;

%% =========================
%  PCA on cleaned data
% ==========================
[coeff, score, latent, tsq, explained] = pca(Xz);  % no NaNs now

%  score   : [Nvox x nRec] PC scores (voxel-wise maps)
%  coeff   : [nRec x nRec] loadings
%  explained: [1 x nRec]   % variance per PC

nPC = size(score, 2);

% Quick sanity check:
fprintf('\nPCA done. First 5 PCs explain (%% variance):\n');
disp(explained(1:min(5,nPC)));

% Optional diagnostic: how many NaNs in scores?
n_nan_scores = sum(isnan(score(:)));
if n_nan_scores > 0
    warning('There are %d NaNs in PC scores even after cleaning.', n_nan_scores);
else
    fprintf('No NaNs in PC scores after cleaning. ✅\n');
end

%% =========================
%  Save PC maps as 1D files
% ==========================
% 1) All PCs in one 1D file (Nvox x nPC)
all_pc_file = fullfile(output_dir, 'PC_all_scores.1D');
dlmwrite(all_pc_file, score, 'delimiter', '\t', 'precision', '%.6f');
fprintf('Saved all PC scores to: %s\n', all_pc_file);

% 2) Each PC in its own 1D file (one column)
for k = 1:nPC
    pc_vec   = score(:,k);              % [Nvox x 1]
    pc_fname = fullfile(output_dir, sprintf('PC%02d.1D', k));
    dlmwrite(pc_fname, pc_vec, 'delimiter', '\t', 'precision', '%.6f');
    fprintf('Saved PC%02d to: %s\n', k, pc_fname);
end

%% =========================
%  Save PCA metadata (loadings etc.)
% ==========================
meta_file = fullfile(output_dir, 'PCA_meta.mat');
save(meta_file, 'coeff', 'latent', 'tsq', 'explained', ...
     'receptor_names', 'mu', 'sigma', 'nPC');
fprintf('Saved PCA metadata to: %s\n', meta_file);

fprintf('\nDone. PC maps + metadata in: %s\n', output_dir);

%% =======================================================================
%                              FUNCTIONS
% ========================================================================

function outdir = ensure_outdir(baseName)
% Make sure we have a real, writable directory
% If a file already exists with that name, append "_dir".

    try
        if exist(baseName, 'file') && ~exist(baseName, 'dir')
            warning('"%s" exists as a file. Using "%s_dir" instead.', baseName, [baseName '_dir']);
            baseName = [baseName '_dir'];
        end
        if ~exist(baseName, 'dir')
            mkdir(baseName);
        end

        % Test write permission
        testFile = fullfile(baseName, sprintf('._writetest_%s.txt', char(java.util.UUID.randomUUID)));
        fid = fopen(testFile, 'w');
        if fid == -1
            error('Cannot write to "%s"', baseName);
        end
        fclose(fid);
        delete(testFile);
    catch
        % Fallback to local folder if scratch dir not writable
        warning('Falling back to local ./output/neuro_PCA_maps directory.');
        baseName = fullfile(pwd, 'output', 'neuro_PCA_maps');
        if ~exist(baseName, 'dir')
            mkdir(baseName);
        end
    end

    outdir = baseName;
end
