%% Pooled SMR-CMC Cross-Correlation Analysis
% This script performs pooled SMR-CMC cross-correlation analysis using
% permutation-significant observations.
%
% Required inputs:
%   smr_timecourse : [360 x 55] SMR time courses
%   cmc_timecourse : [360 x 55] CMC time courses
%
% Data organization:
%   Columns 1:28  = 10% task
%   Columns 29:55 = 15% task
%
% Cross-correlation convention:
%   xcorr(CMC, SMR)
%   Positive lag = SMR leads CMC
%   Negative lag = CMC leads SMR
%
% Outputs are saved to the relative folder:
%   analysis_output/

clearvars -except smr_timecourse cmc_timecourse

assert(exist('smr_timecourse', 'var') == 1, ...
    'Required input "smr_timecourse" was not found.');
assert(exist('cmc_timecourse', 'var') == 1, ...
    'Required input "cmc_timecourse" was not found.');

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

outputDir = fullfile(scriptDir, 'analysis_output');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% Parameters
dt     = 0.05;
t_full = (-8 : dt : 10-dt)';

win_idx = t_full >= 0 & t_full < 7;      % Contraction phase: 0-7 s

MAX_LAG_SEC = 3;
max_lag_pts = round(MAX_LAG_SEC / dt);
lag_sec     = (-max_lag_pts:max_lag_pts)' * dt;
n_lags      = numel(lag_sec);

ALPHA          = 0.05;      % Cluster-level significance threshold
N_PERM         = 1000;      % Number of observation-level permutations
PERM_ALPHA     = 0.05;
N_CLUSTER_PERM = 5000;      % Number of group-level cluster permutations
CLUST_ALPHA    = 0.05;      % Cluster-forming threshold
CLUSTER_SEED   = 42;
PEAK_MODE      = 'positive';
TAIL = 'right';

% Whether negative significant clusters are displayed in the figure.
% When false, negative clusters are not plotted or interpreted.
REPORT_NEG = true;

% The group peak lag is derived from the mean CCF of significant observations.

%% Observation information
ids_10 = setdiff(1:30, [8, 16]);         % N = 28
ids_15 = setdiff(1:30, [3, 5, 16]);      % N = 27

N_10  = numel(ids_10);
N_15  = numel(ids_15);
N_obs = N_10 + N_15;

obs_subject   = [ids_10(:); ids_15(:)];
obs_condition = [repmat({'10%'}, N_10, 1); repmat({'15%'}, N_15, 1)];

fprintf('\nPooled analysis: %d task observations from %d participants\n', ...
    N_obs, numel(unique(obs_subject)));

%% CMC group definition
cmc_pos_10 = [1, 4, 6, 9, 10, 11, 23, 25, 26];
cmc_pos_15 = [6, 9, 10, 11, 23, 25, 26];

is_cmc_pos = false(N_obs, 1);
for k = 1:N_obs
    if strcmp(obs_condition{k}, '10%')
        is_cmc_pos(k) = ismember(obs_subject(k), cmc_pos_10);
    else
        is_cmc_pos(k) = ismember(obs_subject(k), cmc_pos_15);
    end
end
n_pos = sum(is_cmc_pos);
n_neg = sum(~is_cmc_pos);
fprintf('CMC+: %d observations, CMC-: %d observations\n', n_pos, n_neg);

%% Compute observation-level CCF
smr_w = smr_timecourse(win_idx, 1:N_obs);
cmc_w = cmc_timecourse(win_idx, 1:N_obs);

ccf_r = nan(n_lags, N_obs);

for k = 1:N_obs
    smr = smr_w(:, k);
    cmc = cmc_w(:, k);

    if std(smr) == 0 || std(cmc) == 0
        warning('Observation %d has zero variance. Skipped.', k);
        continue
    end

    % xcorr(CMC, SMR): positive lag means SMR leads CMC
    ccf_r(:, k) = xcorr(zscore(cmc), zscore(smr), max_lag_pts, 'coeff');
end

ccf_r_clip = max(min(ccf_r, 0.999999), -0.999999);
ccf_z      = atanh(ccf_r_clip);

%% Observation-level permutation test
individual_peak_r   = nan(N_obs, 1);
individual_peak_lag = nan(N_obs, 1);
individual_perm_p   = nan(N_obs, 1);
individual_sig      = false(N_obs, 1);

rng(2026, 'twister');
fprintf('\nRunning observation-level permutation tests (%d × %d)...\n', N_PERM, N_obs);

for k = 1:N_obs
    this_curve = ccf_r(:, k);

    switch PEAK_MODE
        case 'positive'
            [individual_peak_r(k), idx_k] = max(this_curve);
        case 'abs'
            [~, idx_k] = max(abs(this_curve));
            individual_peak_r(k) = this_curve(idx_k);
    end
    individual_peak_lag(k) = lag_sec(idx_k);

    smr_zk = zscore(smr_w(:, k));
    cmc_zk = zscore(cmc_w(:, k));

    null_peak_k = nan(N_PERM, 1);
    for pi_ = 1:N_PERM
        cmc_perm = cmc_zk(randperm(numel(cmc_zk)));
        r_perm   = xcorr(cmc_perm, smr_zk, max_lag_pts, 'coeff');

        switch PEAK_MODE
            case 'positive', null_peak_k(pi_) = max(r_perm);
            case 'abs',      null_peak_k(pi_) = max(abs(r_perm));
        end
    end

    switch PEAK_MODE
        case 'positive', test_stat = individual_peak_r(k);
        case 'abs',      test_stat = abs(individual_peak_r(k));
    end

    individual_perm_p(k) = (1 + sum(null_peak_k >= test_stat)) / (N_PERM + 1);
    individual_sig(k)    = individual_perm_p(k) < PERM_ALPHA;

    if mod(k, 10) == 0 || k == N_obs
        fprintf('  %d / %d observations done\n', k, N_obs);
    end
end

sig_obs_idx = find(individual_sig);
N_sig_obs   = numel(sig_obs_idx);
fprintf('Permutation-significant observations: %d / %d (p < %.2f)\n', ...
    N_sig_obs, N_obs, PERM_ALPHA);

if N_sig_obs < 3
    error('Fewer than three significant observations are available for group-level analysis.');
end

%% Group mean using permutation-significant observations
ccf_r_sig = ccf_r(:, sig_obs_idx);
ccf_z_sig = ccf_z(:, sig_obs_idx);
subj_sig  = obs_subject(sig_obs_idx);
cmcpos_sig = is_cmc_pos(sig_obs_idx);

n_sig_pos = sum(cmcpos_sig);
n_sig_neg = sum(~cmcpos_sig);
N_sig_subjects = numel(unique(subj_sig));

n_valid_sig = sum(~isnan(ccf_z_sig), 2);
mean_z_sig  = mean(ccf_z_sig, 2, 'omitnan');
se_z_sig    = std(ccf_z_sig, 0, 2, 'omitnan') ./ sqrt(n_valid_sig);
mean_r_sig  = tanh(mean_z_sig);
ci_lo_sig   = tanh(mean_z_sig - tinv(0.975, n_valid_sig - 1) .* se_z_sig);
ci_hi_sig   = tanh(mean_z_sig + tinv(0.975, n_valid_sig - 1) .* se_z_sig);

fprintf('\nGroup-level sample: %d observations (%d subjects) | CMC+: %d, CMC-: %d\n', ...
    N_sig_obs, N_sig_subjects, n_sig_pos, n_sig_neg);

%% Cluster-based permutation test using significant observations
fprintf('\nCluster permutation tail: ''%s''\n', TAIL);

fprintf('Cluster permutation on SIGNIFICANT observations (%d perms)...\n', ...
    N_CLUSTER_PERM);
CL_sig = run_cluster_perm(ccf_z_sig, subj_sig, CLUST_ALPHA, ALPHA, ...
                          N_CLUSTER_PERM, CLUSTER_SEED + 1, TAIL);

print_clusters('Significant observations', CL_sig, lag_sec, ALPHA);

T_cluster_sig = cluster_table(CL_sig, lag_sec, ALPHA);
writetable(T_cluster_sig, fullfile(outputDir, 'SMR_CMC_CCF_cluster_permutation_significant_only.csv'));

%% Peak lag
[peak_r_sig, peak_idx_sig] = pick_peak(mean_r_sig, PEAK_MODE);
peak_lag_sig = lag_sec(peak_idx_sig);

peak_idx = peak_idx_sig;
peak_lag = peak_lag_sig;

% Cluster-level permutation p-value for the cluster containing the peak lag.
[peak_cluster_p_sig, peak_in_cluster_sig] = find_cluster_p_at_lag(CL_sig, peak_idx_sig);

fprintf('\nPeak: r = %.3f at %+.2f s\n', peak_r_sig, peak_lag_sig);
if peak_in_cluster_sig
    fprintf('  Peak is within a significant cluster, cluster p = %.4f\n', peak_cluster_p_sig);
else
    fprintf('  Peak is not within a significant cluster; no point-wise p-value is reported.\n');
end
fprintf('Group peak lag: %+.2f s (derived from significant observations)\n', peak_lag);

r_at_group_peak_lag = ccf_r(peak_idx, :)';

T_obs = table((1:N_obs)', obs_subject, obs_condition, is_cmc_pos, ...
    individual_peak_r, individual_peak_lag, ...
    individual_perm_p, individual_sig, ...
    r_at_group_peak_lag, repmat(peak_lag, N_obs, 1), ...
    'VariableNames', {'Observation', 'SubjectID', 'Condition', 'CMC_positive', ...
                      'IndividualPeak_r', 'IndividualPeakLag_s', ...
                      'Permutation_p', 'Permutation_significant', ...
                      'r_at_GroupPeakLag', 'GroupPeakLag_s'});
writetable(T_obs, fullfile(outputDir, 'SMR_CMC_CCF_observation_level_results.csv'));

%% CMC+ vs CMC- at the group peak lag using a linear mixed-effects model
z_at_peak = ccf_z(peak_idx, :)';
r_at_peak = tanh(z_at_peak);

grp_label = repmat({'CMC-'}, N_obs, 1);
grp_label(is_cmc_pos) = {'CMC+'};

tbl_grp = table(z_at_peak, categorical(obs_subject), ...
    categorical(grp_label, {'CMC-', 'CMC+'}), ...      % CMC- is the reference level
    'VariableNames', {'Z', 'Subject', 'Group'});

lme_grp = fitlme(tbl_grp, 'Z ~ 1 + Group + (1|Subject)');
ct_grp  = lme_grp.Coefficients;

beta_grp = ct_grp.Estimate(2);
se_grp   = ct_grp.SE(2);
t_grp    = ct_grp.tStat(2);
df_grp   = ct_grp.DF(2);
p_grp    = ct_grp.pValue(2);

r_pos = r_at_peak(is_cmc_pos);
r_neg = r_at_peak(~is_cmc_pos);

fprintf('\n===== CMC+ vs CMC- at peak lag (%+.2f s) =====\n', peak_lag);
disp(ct_grp);
fprintf('Delta z = %.3f (SE = %.3f), t(%.1f) = %.3f, p = %.4f\n', ...
    beta_grp, se_grp, df_grp, t_grp, p_grp);
fprintf('  CMC+ (N=%d): M = %.3f, SD = %.3f, Mdn = %.3f\n', ...
    n_pos, mean(r_pos,'omitnan'), std(r_pos,'omitnan'), median(r_pos,'omitnan'));
fprintf('  CMC- (N=%d): M = %.3f, SD = %.3f, Mdn = %.3f\n', ...
    n_neg, mean(r_neg,'omitnan'), std(r_neg,'omitnan'), median(r_neg,'omitnan'));

%% Save numeric results
T_lag = table(lag_sec, mean_r_sig, ci_lo_sig, ci_hi_sig, ...
    CL_sig.t_obs, CL_sig.sig, CL_sig.sig_pos, CL_sig.sig_neg, ...
    'VariableNames', {'Lag_s', ...
        'Mean_r_sig','CI_low_sig','CI_high_sig','Cluster_t_sig','Cluster_sig_sig', ...
        'Cluster_sig_sig_positive','Cluster_sig_sig_negative'});
writetable(T_lag, fullfile(outputDir, 'SMR_CMC_CCF_pooled_lag_results.csv'));

T_peak_group = table({'CMC-'; 'CMC+'}, [n_neg; n_pos], ...
    [mean(r_neg,'omitnan'); mean(r_pos,'omitnan')], ...
    [std(r_neg,'omitnan');  std(r_pos,'omitnan')], ...
    [median(r_neg,'omitnan'); median(r_pos,'omitnan')], ...
    [NaN; beta_grp], [NaN; t_grp], [NaN; df_grp], [NaN; p_grp], ...
    'VariableNames', {'Group','N','Mean_r','SD_r','Median_r', ...
                      'LMM_beta_z','LMM_t','LMM_df','LMM_p'});
writetable(T_peak_group, fullfile(outputDir, 'SMR_CMC_CCF_CMCgroup_at_peak_lag.csv'));

save(fullfile(outputDir, 'SMR_CMC_CCF_pooled_results.mat'), ...
    'ccf_r', 'ccf_z', 'lag_sec', ...
    'ccf_r_sig', 'ccf_z_sig', 'sig_obs_idx', 'subj_sig', 'cmcpos_sig', ...
    'mean_r_sig', 'ci_lo_sig', 'ci_hi_sig', ...
    'N_sig_obs', 'N_sig_subjects', 'n_sig_pos', 'n_sig_neg', ...
    'CL_sig', 'T_cluster_sig', ...
    'N_CLUSTER_PERM', 'CLUST_ALPHA', 'ALPHA', 'TAIL', 'REPORT_NEG', ...
    'peak_r_sig', 'peak_lag_sig', 'peak_idx_sig', 'peak_lag', 'peak_idx', ...
    'peak_cluster_p_sig', 'peak_in_cluster_sig', ...
    'obs_subject', 'obs_condition', 'is_cmc_pos', 'r_at_group_peak_lag', ...
    'individual_peak_r', 'individual_peak_lag', 'individual_perm_p', ...
    'individual_sig', 'T_obs', ...
    'lme_grp', 'beta_grp', 'se_grp', 't_grp', 'df_grp', 'p_grp', ...
    'N_PERM', 'PERM_ALPHA');

%% Summary figure
C_POS      = [0.22 0.51 0.82];      % CMC+
C_NEG      = [0.72 0.72 0.72];      % CMC-
C_IND      = [0.80 0.80 0.80];      % individual curves
C_CLUS_POS = [0.90 0.30 0.30];      % Positive significant cluster
C_CLUS_NEG = [0.30 0.45 0.78];      % Negative significant cluster

fig = figure('Color', 'w', 'Position', [80 80 1180 420]);

% Left panel: CMC- versus CMC+ at the group peak lag using all observations.
subplot(1, 3, 1); hold on;

rng(42, 'twister');
jit = 0.22;

draw_box(1, r_neg, 0.55, C_NEG, [0 0 0]);
draw_box(2, r_pos, 0.55, C_POS, [0 0 0]);

scatter(1 + (rand(n_neg,1)-.5)*jit, r_neg, 22, 'k', 'filled', 'MarkerFaceAlpha', 0.85);
scatter(2 + (rand(n_pos,1)-.5)*jit, r_pos, 22, 'k', 'filled', 'MarkerFaceAlpha', 0.85);

y_max = max([r_neg; r_pos], [], 'omitnan') + 0.14;
plot([1 2], [y_max y_max], 'k-', 'LineWidth', 1.2);

if     p_grp < 0.001, sig_str = '***';
elseif p_grp < 0.01,  sig_str = '**';
elseif p_grp < 0.05,  sig_str = '*';
else,                 sig_str = sprintf('p = %.3f', p_grp);
end
text(1.5, y_max + 0.02, sig_str, 'HorizontalAlignment', 'center', ...
    'FontSize', 16, 'FontWeight', 'bold');

xlim([0.4 2.6]); ylim([-1.05 1.05]);
xticks([1 2]); xticklabels({'CMC-', 'CMC+'});
yticks(-1:0.5:1);
ylabel('Cross-correlation coefficient', 'FontSize', 12);
set(gca, 'FontSize', 12, 'LineWidth', 1.0, 'TickDir', 'out', 'Box', 'off');

% Right panel: pooled CCF using permutation-significant observations only.
subplot(1, 3, [2 3]); hold on;

for jj = 1:N_sig_obs
    plot(lag_sec, ccf_r_sig(:, jj), '-', 'Color', C_IND, 'LineWidth', 0.5, ...
        'HandleVisibility', 'off');
end

% Mean 95% confidence interval across the full lag range.
h_ci = fill([lag_sec; flipud(lag_sec)], [ci_hi_sig; flipud(ci_lo_sig)], ...
    [0.5 0.5 0.5], 'FaceAlpha', 0.20, 'EdgeColor', 'none');

% Significant clusters.
h_clu_pos = gobjects(0);
h_clu_neg = gobjects(0);

segs_pos = find_segments(CL_sig.sig_pos);
for ss = 1:size(segs_pos, 1)
    ii = segs_pos(ss,1):segs_pos(ss,2);
    hf = fill([lag_sec(ii); flipud(lag_sec(ii))], [ci_hi_sig(ii); flipud(ci_lo_sig(ii))], ...
        C_CLUS_POS, 'FaceAlpha', 0.35, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    if ss == 1, h_clu_pos = hf; end
end

if REPORT_NEG
    segs_neg = find_segments(CL_sig.sig_neg);
    for ss = 1:size(segs_neg, 1)
        ii = segs_neg(ss,1):segs_neg(ss,2);
        hf = fill([lag_sec(ii); flipud(lag_sec(ii))], [ci_hi_sig(ii); flipud(ci_lo_sig(ii))], ...
            C_CLUS_NEG, 'FaceAlpha', 0.30, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        if ss == 1, h_clu_neg = hf; end
    end
end

h_mean = plot(lag_sec, mean_r_sig, 'k-', 'LineWidth', 2.6);

h_pos = gobjects(0);
h_neg = gobjects(0);
if n_sig_pos > 0
    h_pos = plot(lag_sec, tanh(mean(ccf_z_sig(:,  cmcpos_sig), 2, 'omitnan')), '--', ...
        'Color', C_POS, 'LineWidth', 2.0);
end
if n_sig_neg > 0
    h_neg = plot(lag_sec, tanh(mean(ccf_z_sig(:, ~cmcpos_sig), 2, 'omitnan')), '--', ...
        'Color', [0.45 0.45 0.45], 'LineWidth', 2.0);
end
h_ind = plot(nan, nan, '-', 'Color', C_IND, 'LineWidth', 0.8);

plot(peak_lag_sig, peak_r_sig, 'o', 'MarkerFaceColor', [1 0 0], ...
    'MarkerEdgeColor', [1 0 0], 'MarkerSize', 8, 'HandleVisibility', 'off');

if peak_in_cluster_sig
    text(peak_lag_sig + 0.12, peak_r_sig + 0.05, sprintf('cluster p = %.3f', peak_cluster_p_sig), ...
        'FontSize', 9, 'Color', [0.1 0.1 0.1]);
end

xline(0, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0, 'HandleVisibility', 'off');
yline(0, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8, 'HandleVisibility', 'off');

leg_h = [h_mean h_ci];  leg_s = {'Mean', 'Mean 95% CI'};
if ~isempty(h_pos), leg_h(end+1) = h_pos; leg_s{end+1} = sprintf('CMC+ mean (n=%d)', n_sig_pos); end
if ~isempty(h_neg), leg_h(end+1) = h_neg; leg_s{end+1} = sprintf('CMC- mean (n=%d)', n_sig_neg); end
leg_h(end+1) = h_ind;  leg_s{end+1} = 'Individuals';
if ~isempty(h_clu_pos), leg_h(end+1) = h_clu_pos; leg_s{end+1} = 'Positive cluster p < .05'; end
if ~isempty(h_clu_neg), leg_h(end+1) = h_clu_neg; leg_s{end+1} = 'Negative cluster p < .05'; end

xlim([-MAX_LAG_SEC MAX_LAG_SEC]);
xticks(-3:1:3);
xlabel({'Lag (s)', '\leftarrowCMC leads | SMR leads\rightarrow'}, 'FontSize', 12);
ylabel('Cross-correlation coefficient', 'FontSize', 12);
legend(leg_h, leg_s, 'Location', 'northeast', 'FontSize', 9, 'Box', 'off');
set(gca, 'FontSize', 12, 'LineWidth', 1.0, 'TickDir', 'out', 'Box', 'off');

print(fig, fullfile(outputDir, 'SMR_CMC_CCF_summary'), '-dpng', '-r300');
savefig(fig, fullfile(outputDir, 'SMR_CMC_CCF_summary.fig'));

%% Saved files
fprintf('\nSaved files:\n');
fprintf('  SMR_CMC_CCF_pooled_lag_results.csv\n');
fprintf('  SMR_CMC_CCF_cluster_permutation_significant_only.csv\n');
fprintf('  SMR_CMC_CCF_observation_level_results.csv\n');
fprintf('  SMR_CMC_CCF_CMCgroup_at_peak_lag.csv\n');
fprintf('  SMR_CMC_CCF_pooled_results.mat\n');
fprintf('  SMR_CMC_CCF_summary.png / .fig\n');

%% Local functions
function S = run_cluster_perm(z_mat, subj_ids, clust_alpha, alpha, n_perm, seed, tail)
% Cluster permutation: average within subject, one-sample t statistic, and sign-flip permutation.
%   tail = 'two'   : threshold = tinv(1-alpha/2); test positive and negative clusters; null = max|mass|.
%   tail = 'right' : threshold = tinv(1-alpha); test positive clusters only; null = max(mass).
    if nargin < 7, tail = 'right'; end
    n_lags = size(z_mat, 1);
    u = unique(subj_ids);
    N = numel(u);

    S.subject_ids = u;
    S.N_subjects  = N;
    S.subject_z   = nan(n_lags, N);
    for si = 1:N
        S.subject_z(:, si) = mean(z_mat(:, subj_ids == u(si)), 2, 'omitnan');
    end

    S.t_obs   = nan(n_lags, 1);
    S.sig     = false(n_lags, 1);
    S.sig_pos = false(n_lags, 1);
    S.sig_neg = false(n_lags, 1);
    S.mass = []; S.segs = []; S.sign = []; S.p = []; S.null_max = [];

    S.tail = tail;

    if N < 3
        warning('Only %d subjects are available; cluster permutation is skipped.', N);
        S.tcrit = NaN;
        return
    end

    S.t_obs = mean(S.subject_z, 2) ./ (std(S.subject_z, 0, 2) / sqrt(N));

    switch lower(tail)
        case 'two'
            S.tcrit  = tinv(1 - clust_alpha/2, N - 1);
            find_fun = @(t) cluster_mass_two_tailed(t, S.tcrit);
            null_fun = @(m) max(abs(m));
        case 'right'
            S.tcrit  = tinv(1 - clust_alpha, N - 1);
            find_fun = @(t) cluster_mass_right_only(t, S.tcrit);
            null_fun = @(m) max(m);
        otherwise
            error('TAIL must be ''two'' or ''right''.');
    end

    [S.mass, S.segs, S.sign] = find_fun(S.t_obs);

    rng(seed, 'twister');
    S.null_max = zeros(n_perm, 1);
    for ip = 1:n_perm
        flip   = 2 * (rand(1, N) > 0.5) - 1;
        z_perm = S.subject_z .* flip;
        t_perm = mean(z_perm, 2) ./ (std(z_perm, 0, 2) / sqrt(N));
        m = find_fun(t_perm);
        if ~isempty(m), S.null_max(ip) = null_fun(m); end
    end

    S.p = nan(numel(S.mass), 1);
    for ci = 1:numel(S.mass)
        S.p(ci) = (1 + sum(S.null_max >= abs(S.mass(ci)))) / (n_perm + 1);
        if S.p(ci) < alpha
            idx = S.segs(ci,1):S.segs(ci,2);
            S.sig(idx) = true;
            if S.sign(ci) > 0, S.sig_pos(idx) = true; else, S.sig_neg(idx) = true; end
        end
    end
end

function [masses, segs, signs_out] = cluster_mass_two_tailed(t_vec, tcrit)
% Two-tailed supra-threshold clusters; cluster mass is the signed sum of t values.
    masses = []; segs = []; signs_out = [];

    for s = [1, -1]
        idx = find(s * t_vec > tcrit);
        if isempty(idx), continue; end

        bp    = find(diff(idx) > 1);
        seg_s = [idx(1); idx(bp + 1)];
        seg_e = [idx(bp); idx(end)];
        m     = arrayfun(@(a,b) sum(t_vec(a:b)), seg_s, seg_e);

        masses    = [masses; m(:)];                       %#ok<AGROW>
        segs      = [segs; [seg_s(:), seg_e(:)]];         %#ok<AGROW>
        signs_out = [signs_out; repmat(s, numel(m), 1)];  %#ok<AGROW>
    end

    if ~isempty(segs)
        [~, ord] = sort(segs(:,1));
        masses    = masses(ord);
        segs      = segs(ord, :);
        signs_out = signs_out(ord);
    end
end

function [p_at_lag, in_cluster] = find_cluster_p_at_lag(CL, idx)
% Return the p-value of the significant cluster containing the specified lag index.
    p_at_lag   = NaN;
    in_cluster = false;
    if isempty(CL.mass), return; end
    for ci = 1:numel(CL.mass)
        if idx >= CL.segs(ci,1) && idx <= CL.segs(ci,2) && CL.p(ci) < 0.05
            p_at_lag   = CL.p(ci);
            in_cluster = true;
            return
        end
    end
end

function [masses, segs, signs_out] = cluster_mass_right_only(t_vec, tcrit)
% Right-tailed clustering: consecutive samples above +tcrit form a cluster.
    idx = find(t_vec > tcrit);
    if isempty(idx)
        masses = []; segs = []; signs_out = [];
        return
    end

    bp    = find(diff(idx) > 1);
    seg_s = [idx(1); idx(bp + 1)];
    seg_e = [idx(bp); idx(end)];

    masses    = arrayfun(@(a,b) sum(t_vec(a:b)), seg_s, seg_e);
    masses    = masses(:);
    segs      = [seg_s(:), seg_e(:)];
    signs_out = ones(numel(masses), 1);
end

function T = cluster_table(S, lag_sec, alpha)
    if isempty(S.mass)
        T = table();
        return
    end
    dir_col = repmat({'positive'}, numel(S.mass), 1);
    dir_col(S.sign < 0) = {'negative'};
    T = table(lag_sec(S.segs(:,1)), lag_sec(S.segs(:,2)), dir_col, ...
        S.mass(:), S.p(:), S.p(:) < alpha, ...
        'VariableNames', {'Lag_start_s','Lag_end_s','Direction', ...
                          'ClusterMass','p_perm','Significant'});
end

function print_clusters(label, S, lag_sec, alpha)
    fprintf('\n--- Clusters: %s (N = %d subjects) ---\n', label, S.N_subjects);
    if isempty(S.mass)
        fprintf('  No lag samples exceeded the cluster-forming threshold.\n');
        return
    end
    for ci = 1:numel(S.mass)
        if S.sign(ci) > 0, d = 'positive'; else, d = 'negative'; end
        if S.p(ci) < alpha, mk = '  *'; else, mk = ''; end
        fprintf('  %+.2f to %+.2f s [%s]: mass = %+.1f, p = %.4f%s\n', ...
            lag_sec(S.segs(ci,1)), lag_sec(S.segs(ci,2)), d, S.mass(ci), S.p(ci), mk);
    end
    if ~any(S.p < alpha)
        fprintf('  No cluster reached p < %.2f.\n', alpha);
    end
end

function [pk, idx] = pick_peak(curve, mode)
    switch mode
        case 'positive'
            [pk, idx] = max(curve);
        case 'abs'
            [~, idx] = max(abs(curve));
            pk = curve(idx);
        otherwise
            error('PEAK_MODE must be ''positive'' or ''abs''.');
    end
end

function segs = find_segments(mask)
% Return [start, end] indices for consecutive true segments.
    mask = logical(mask(:));
    idx  = find(mask);
    if isempty(idx), segs = []; return; end

    bp    = find(diff(idx) > 1);
    seg_s = [idx(1); idx(bp + 1)];
    seg_e = [idx(bp); idx(end)];
    segs  = [seg_s, seg_e];
end

function draw_box(x, y, w, faceCol, edgeCol)
% Draw a Tukey-style box plot: Q1-Q3 box and whiskers within 1.5 x IQR.
    y = y(~isnan(y));
    if isempty(y), return; end

    q1 = prctile(y, 25);
    q2 = median(y);
    q3 = prctile(y, 75);
    iqr_val = q3 - q1;

    lo_w = min(y(y >= q1 - 1.5*iqr_val));
    hi_w = max(y(y <= q3 + 1.5*iqr_val));

    fill([x-w/2 x+w/2 x+w/2 x-w/2], [q1 q1 q3 q3], faceCol, ...
        'EdgeColor', edgeCol, 'LineWidth', 1.1, 'HandleVisibility', 'off');
    plot([x-w/2 x+w/2], [q2 q2], '-', 'Color', edgeCol, 'LineWidth', 1.8, ...
        'HandleVisibility', 'off');
    plot([x x], [q3 hi_w], '-', 'Color', edgeCol, 'LineWidth', 1.1, ...
        'HandleVisibility', 'off');
    plot([x x], [q1 lo_w], '-', 'Color', edgeCol, 'LineWidth', 1.1, ...
        'HandleVisibility', 'off');
    plot([x-w/4 x+w/4], [hi_w hi_w], '-', 'Color', edgeCol, 'LineWidth', 1.1, ...
        'HandleVisibility', 'off');
    plot([x-w/4 x+w/4], [lo_w lo_w], '-', 'Color', edgeCol, 'LineWidth', 1.1, ...
        'HandleVisibility', 'off');
end