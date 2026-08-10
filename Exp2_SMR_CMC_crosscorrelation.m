%% SMR-CMC Cross-Correlation Across Tasks
% GitHub-ready analysis script
%
% Description:
% This script compares SMR-CMC cross-correlation functions across four
% within-participant tasks: St, 1sin, 3sin, and 4sin.
%
% Required workspace variables:
%   smr_1, cmc_1
%   smr_2, cmc_2
%   smr_3, cmc_3
%   smr_4, cmc_4
%
% Each matrix must have dimensions [360 x N], where columns correspond
% to participants. The same participant must use the same subject ID
% across all four tasks.
%
% Cross-correlation convention:
%   xcorr(CMC, SMR)
%   Positive lag: SMR leads CMC
%   Negative lag: CMC leads SMR
%
% Main analyses:
%   1. Compute observation-level SMR-CMC cross-correlation functions.
%   2. Test positive lag-wise correlations using a one-tailed
%      cluster-based permutation test.
%   3. Determine a task-specific peak lag from each group-mean CCF.
%   4. Compare individual r values at those task-specific peak lags
%      using a Friedman test with Dunn's post-hoc comparisons.
%   5. Examine individual peak-lag distributions as a timing diagnostic.
%
clearvars -except smr_1 smr_2 smr_3 smr_4 cmc_1 cmc_2 cmc_3 cmc_4

%% User configuration
SMR_ALL = {smr_1, smr_2, smr_3, smr_4};
CMC_ALL = {cmc_1, cmc_2, cmc_3, cmc_4};

TASK_NAMES = {'St', '1sin', '3sin', '4sin'};

% The same participant must use the same ID across all four tasks.
SUBJ_IDS = { 1:9, 1:9, 1:9, 1:9 };

dt          = 0.05;
t_full      = (-8 : dt : 10-dt)';
win_idx     = t_full >= 0 & t_full < 7;    % contraction phase: 0–7 s
MAX_LAG_SEC = 3;                            % Lag range: +/-3 s
ALPHA       = 0.05;
PEAK_MODE   = 'positive';                   % Peak definition for each task

% Cluster-based permutation parameters
N_PERM      = 5000;
CLUST_ALPHA = 0.05;      % Pointwise cluster-forming threshold, one-tailed
PERM_SEED   = 42;

% Plot colors
C_TASK = [0.851 0.851 0.851;    % St    #D9D9D9
          0.886 0.353 0.290;    % 1sin  #E25A4A
          0.427 0.784 0.604;    % 3sin  #6DC89A
          0.318 0.549 0.612];   % 4sin  #518C9C
C_LINE = C_TASK;
C_LINE(1, :) = [0.45 0.45 0.45];   % Darker gray for line visibility
C_PEAK = [0.85 0.10 0.10];         % Peak marker color

OUT_PREFIX = 'SMR_CMC_CCF_taskcomp';

% Save all outputs relative to this script.
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

OUTPUT_DIR = fullfile(scriptDir, 'analysis_output');
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

%% Derived parameters and input validation
n_task      = numel(SMR_ALL);
max_lag_pts = round(MAX_LAG_SEC / dt);
lag_sec     = (-max_lag_pts:max_lag_pts)' * dt;
n_lags      = numel(lag_sec);

assert(numel(CMC_ALL) == n_task && numel(TASK_NAMES) == n_task && ...
       numel(SUBJ_IDS) == n_task, 'The number of task inputs is inconsistent.');

for tk = 1:n_task
    assert(size(SMR_ALL{tk}, 2) == size(CMC_ALL{tk}, 2), ...
        '%s: SMR and CMC must have the same number of columns.', TASK_NAMES{tk});
    assert(size(SMR_ALL{tk}, 2) == numel(SUBJ_IDS{tk}), ...
        '%s: SUBJ_IDS length must match the number of data columns.', TASK_NAMES{tk});
    assert(size(SMR_ALL{tk}, 1) == numel(t_full), ...
        '%s: the number of rows must match t_full.', TASK_NAMES{tk});
end

%% Compute observation-level cross-correlation functions
ccf_r_all   = [];
obs_subject = [];
obs_task    = {};

for tk = 1:n_task

    smr_w = SMR_ALL{tk}(win_idx, :);
    cmc_w = CMC_ALL{tk}(win_idx, :);
    n_col = size(smr_w, 2);

    ccf_this = nan(n_lags, n_col);

    for k = 1:n_col
        smr = smr_w(:, k);
        cmc = cmc_w(:, k);

        ok  = ~isnan(smr) & ~isnan(cmc);
        smr = smr(ok);
        cmc = cmc(ok);

        if numel(smr) <= max_lag_pts + 1
            warning('%s observation %d: too few valid samples; skipped.', TASK_NAMES{tk}, k);
            continue
        end
        if std(smr) == 0 || std(cmc) == 0
            warning('%s observation %d: zero variance; skipped.', TASK_NAMES{tk}, k);
            continue
        end

        % xcorr(CMC, SMR): positive lag means that SMR leads CMC.
        r = xcorr(zscore(cmc), zscore(smr), max_lag_pts, 'coeff');
        ccf_this(:, k) = r(:);
    end

    ccf_r_all   = [ccf_r_all, ccf_this];                              %#ok<AGROW>
    obs_subject = [obs_subject; SUBJ_IDS{tk}(:)];                     %#ok<AGROW>
    obs_task    = [obs_task; repmat(TASK_NAMES(tk), n_col, 1)];       %#ok<AGROW>
end

N_obs = size(ccf_r_all, 2);

ccf_r_clip = max(min(ccf_r_all, 0.999999), -0.999999);
ccf_z      = atanh(ccf_r_clip);          % Fisher z transform for cluster inference

subj_u = unique(obs_subject);
N_subj = numel(subj_u);

fprintf('\nTask comparison: %d observations, %d participants, %d tasks\n', ...
    N_obs, N_subj, n_task);
for tk = 1:n_task
    fprintf('  %s: N = %d\n', TASK_NAMES{tk}, sum(strcmp(obs_task, TASK_NAMES{tk})));
end

%% Build participant-by-task column index table
col_idx = nan(N_subj, n_task);
for si = 1:N_subj
    for tk = 1:n_task
        idx = find(obs_subject == subj_u(si) & strcmp(obs_task, TASK_NAMES{tk}), 1);
        if ~isempty(idx), col_idx(si, tk) = idx; end
    end
end

complete_mask = all(~isnan(col_idx), 2);
for si = 1:N_subj
    if ~complete_mask(si), continue; end
    if any(any(isnan(ccf_r_all(:, col_idx(si, :)))))
        complete_mask(si) = false;
    end
end

subj_c    = subj_u(complete_mask);
col_idx_c = col_idx(complete_mask, :);
n_rm      = numel(subj_c);

if sum(~complete_mask) > 0
    fprintf('\n[Friedman] %d participants had incomplete data and were excluded (IDs: %s)\n', ...
        sum(~complete_mask), num2str(subj_u(~complete_mask)'));
end
fprintf('[Friedman] Complete data: N = %d participants x %d tasks\n', n_rm, n_task);
assert(n_rm >= 3, 'Too few participants with complete data.');

%% Lag-wise group means, confidence intervals, and one-tailed t statistics
task_mean_z = nan(n_lags, n_task);
task_se_z   = nan(n_lags, n_task);
task_t0     = nan(n_lags, n_task);
task_p0     = nan(n_lags, n_task);    % Uncorrected one-tailed p-values for descriptive output
task_n      = nan(n_task, 1);
ci_lo_r     = nan(n_lags, n_task);
ci_hi_r     = nan(n_lags, n_task);

for tk = 1:n_task
    sel  = strcmp(obs_task, TASK_NAMES{tk});
    Zt   = ccf_z(:, sel);
    Zt   = Zt(:, all(~isnan(Zt), 1));
    n_tk = size(Zt, 2);
    task_n(tk) = n_tk;
    if n_tk < 3, continue; end

    task_mean_z(:, tk) = mean(Zt, 2);
    task_se_z(:, tk)   = std(Zt, 0, 2) / sqrt(n_tk);
    task_t0(:, tk)     = task_mean_z(:, tk) ./ task_se_z(:, tk);
    task_p0(:, tk)     = 1 - tcdf(task_t0(:, tk), n_tk - 1);

    tcrit = tinv(1 - ALPHA/2, n_tk - 1);
    ci_lo_r(:, tk) = tanh(task_mean_z(:, tk) - tcrit * task_se_z(:, tk));
    ci_hi_r(:, tk) = tanh(task_mean_z(:, tk) + tcrit * task_se_z(:, tk));
end

task_mean_r = tanh(task_mean_z);

%% Determine the peak lag of the group-mean CCF for each task
task_peak_idx = nan(n_task, 1);
task_peak_lag = nan(n_task, 1);
task_peak_r   = nan(n_task, 1);

for tk = 1:n_task
    curve = task_mean_r(:, tk);
    switch PEAK_MODE
        case 'positive'
            [task_peak_r(tk), idx_t] = max(curve);
        case 'abs'
            [~, idx_t] = max(abs(curve));
            task_peak_r(tk) = curve(idx_t);
        otherwise
            error('PEAK_MODE must be ''positive'' or ''abs''.');
    end
    task_peak_idx(tk) = idx_t;
    task_peak_lag(tk) = lag_sec(idx_t);
end

%% Cluster-based permutation test for r > 0 across lags
rng(PERM_SEED);

task_sig0  = false(n_lags, n_task);
clust_tbls = cell(n_task, 1);

fprintf('\n===== Lag-wise r > 0: one-tailed cluster permutation (%d permutations) =====\n', N_PERM);

for tk = 1:n_task

    sel  = strcmp(obs_task, TASK_NAMES{tk});
    Zt   = ccf_z(:, sel);
    Zt   = Zt(:, all(~isnan(Zt), 1));
    n_tk = size(Zt, 2);
    if n_tk < 3
        fprintf('%s: too few valid participants; skipped.\n', TASK_NAMES{tk});
        clust_tbls{tk} = table();
        continue
    end

    tcrit_c = tinv(1 - CLUST_ALPHA, n_tk - 1);

    t_obs = mean(Zt, 2) ./ (std(Zt, 0, 2) / sqrt(n_tk));
    [obs_mass, obs_segs] = cluster_mass_right(t_obs, tcrit_c);

    max_mass = zeros(N_PERM, 1);
    for ip = 1:N_PERM
        sgn = sign(rand(1, n_tk) - 0.5);
        Zp  = Zt .* sgn;
        t_p = mean(Zp, 2) ./ (std(Zp, 0, 2) / sqrt(n_tk));
        mp  = cluster_mass_right(t_p, tcrit_c);
        if ~isempty(mp), max_mass(ip) = max(mp); end
    end

    fprintf('\n%s (N = %d, pointwise threshold t > %.3f):\n', TASK_NAMES{tk}, n_tk, tcrit_c);
    if isempty(obs_mass)
        fprintf('  No supra-threshold clusters.\n');
        clust_tbls{tk} = table();
    else
        p_clust = arrayfun(@(m) (1 + sum(max_mass >= m)) / (N_PERM + 1), obs_mass);
        for c = 1:numel(obs_mass)
            i1 = obs_segs(c,1); i2 = obs_segs(c,2);
            [tmax, itmax] = max(t_obs(i1:i2));
            fprintf('  %+.2f s ~ %+.2f s | mass = %.1f, t_max = %.2f @ %+.2f s, p = %.4f%s\n', ...
                lag_sec(i1), lag_sec(i2), obs_mass(c), tmax, ...
                lag_sec(i1 + itmax - 1), p_clust(c), star_str(p_clust(c)));
            if p_clust(c) < ALPHA
                task_sig0(i1:i2, tk) = true;
            end
        end
        clust_tbls{tk} = table( ...
            repmat(TASK_NAMES(tk), numel(obs_mass), 1), ...
            lag_sec(obs_segs(:,1)), lag_sec(obs_segs(:,2)), ...
            obs_mass(:), p_clust(:), p_clust(:) < ALPHA, ...
            'VariableNames', {'Task', 'Lag_start_s', 'Lag_end_s', ...
                              'ClusterMass', 'p_perm', 'Significant'});
    end
end

T_clust = vertcat(clust_tbls{:});
if ~isempty(T_clust)
    writetable(T_clust, fullfile(OUTPUT_DIR, [OUT_PREFIX '_cluster_permutation.csv']));
end

%% Compare tasks using individual r values at task-specific peak lags
r_at_taskpeak = nan(N_obs, 1);
for k = 1:N_obs
    tk = find(strcmp(TASK_NAMES, obs_task{k}), 1);
    r_at_taskpeak(k) = ccf_r_all(task_peak_idx(tk), k);
end

% Wide-format matrix for participants with complete data
Rw = nan(n_rm, n_task);
for si = 1:n_rm
    Rw(si, :) = r_at_taskpeak(col_idx_c(si, :))';
end

fprintf('\n===== Task-specific peak lags and individual r values =====\n');
for tk = 1:n_task
    fprintf('  %-6s peak lag = %+.2f s (group r = %.3f) | Individual r: median = %.3f, IQR = %.3f–%.3f, >0: %d/%d\n', ...
        TASK_NAMES{tk}, task_peak_lag(tk), task_peak_r(tk), ...
        median(Rw(:, tk)), prctile(Rw(:, tk), 25), prctile(Rw(:, tk), 75), ...
        sum(Rw(:, tk) > 0), n_rm);
end

%% Friedman test and Dunn's post-hoc comparisons
FD = friedman_dunn(Rw, TASK_NAMES);

fprintf('\n===== Friedman test: r at task-specific peak lag =====\n');
fprintf('N = %d participants, k = %d tasks\n', FD.n, FD.k);
fprintf('Friedman chi2(%d) = %.3f, p = %.4f%s\n', ...
    FD.df, FD.chi2, FD.p, star_str(FD.p));
fprintf("Kendall's W = %.3f\n", FD.W);
fprintf('Task rank sums: ');
for tk = 1:n_task
    fprintf('%s = %.1f  ', TASK_NAMES{tk}, FD.rank_sum(tk));
end
fprintf('\n');

fprintf("\nDunn's multiple comparisons (Bonferroni correction, %d comparisons):\n", FD.n_pair);
fprintf('%-6s %-6s %10s %8s %12s %12s\n', 'A', 'B', 'RankSumDiff', 'z', 'p_unadj', 'p_adj');
for pp = 1:FD.n_pair
    fprintf('%-6s %-6s %10.2f %8.3f %12.4f %12.4f%s\n', ...
        FD.pair_A{pp}, FD.pair_B{pp}, FD.rank_diff(pp), FD.z(pp), ...
        FD.p_unadj(pp), FD.p_adj(pp), star_str(FD.p_adj(pp)));
end

pair_list = FD.pair_list;
n_pair    = FD.n_pair;
p_dunn    = FD.p_adj;

writetable(FD.T_pair, fullfile(OUTPUT_DIR, [OUT_PREFIX '_dunn_at_task_peak_lag.csv']));

%% Individual peak-lag distribution for timing-consistency diagnostics
ind_peak_lag = nan(N_obs, 1);
for k = 1:N_obs
    curve = ccf_r_all(:, k);
    if all(isnan(curve)), continue; end
    switch PEAK_MODE
        case 'positive', [~, idx_k] = max(curve);
        case 'abs',      [~, idx_k] = max(abs(curve));
    end
    ind_peak_lag(k) = lag_sec(idx_k);
end

Lw = nan(n_rm, n_task);
for si = 1:n_rm
    Lw(si, :) = ind_peak_lag(col_idx_c(si, :))';
end

FD_lag = friedman_dunn(Lw, TASK_NAMES);
fprintf('\n[Diagnostic] Individual peak-lag task difference: chi2(%d) = %.3f, p = %.4f, W = %.3f\n', ...
    FD_lag.df, FD_lag.chi2, FD_lag.p, FD_lag.W);
fprintf('        SD of individual peak lag by task: ');
sd_lag = std(Lw, 0, 1);
for tk = 1:n_task
    fprintf('%s = %.2f s  ', TASK_NAMES{tk}, sd_lag(tk));
end
fprintf('\n');

%% Save numeric results
T_lag = table(lag_sec, 'VariableNames', {'Lag_s'});
for tk = 1:n_task
    nm = matlab.lang.makeValidName(TASK_NAMES{tk});
    T_lag.(['Mean_r_' nm])   = task_mean_r(:, tk);
    T_lag.(['CI_lo_' nm])    = ci_lo_r(:, tk);
    T_lag.(['CI_hi_' nm])    = ci_hi_r(:, tk);
    T_lag.(['t_' nm])        = task_t0(:, tk);
    T_lag.(['p_1tail_' nm])  = task_p0(:, tk);
    T_lag.(['SigClust_' nm]) = task_sig0(:, tk);
end
writetable(T_lag, fullfile(OUTPUT_DIR, [OUT_PREFIX '_lagwise.csv']));

obs_task_peak_lag = nan(N_obs, 1);
for k = 1:N_obs
    obs_task_peak_lag(k) = task_peak_lag(strcmp(TASK_NAMES, obs_task{k}));
end

T_obs = table((1:N_obs)', obs_subject, obs_task, obs_task_peak_lag, ...
    r_at_taskpeak, ind_peak_lag, ...
    'VariableNames', {'Observation', 'SubjectID', 'Task', 'TaskPeakLag_s', ...
                      'r_at_TaskPeakLag', 'IndividualPeakLag_s'});
writetable(T_obs, fullfile(OUTPUT_DIR, [OUT_PREFIX '_observation_level.csv']));

save(fullfile(OUTPUT_DIR, [OUT_PREFIX '_results.mat']), ...
    'ccf_r_all', 'ccf_z', 'lag_sec', 'obs_subject', 'obs_task', 'TASK_NAMES', ...
    'task_mean_r', 'task_se_z', 'task_t0', 'task_p0', 'ci_lo_r', 'ci_hi_r', ...
    'task_sig0', 'clust_tbls', 'T_clust', ...
    'task_peak_idx', 'task_peak_lag', 'task_peak_r', ...
    'r_at_taskpeak', 'Rw', 'Lw', 'subj_c', 'FD', 'FD_lag', 'T_obs');

%% Figure 1: task CCF curves with significant clusters emphasized
LW_NS  = 1.3;
LW_SIG = 3.4;
A_NS   = 0.45;

fig1 = figure('Color', 'w', 'Position', [80 80 840 560]); hold on;

for tk = 1:n_task
    fill([lag_sec; flipud(lag_sec)], [ci_hi_r(:, tk); flipud(ci_lo_r(:, tk))], ...
        C_LINE(tk, :), 'FaceAlpha', 0.10, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
end

xline(0, 'k--', 'LineWidth', 1.0, 'Alpha', 0.4, 'HandleVisibility', 'off');
yline(0, 'k-',  'LineWidth', 0.6, 'Alpha', 0.25, 'HandleVisibility', 'off');

h = gobjects(1, n_task);
for tk = 1:n_task
    h(tk) = plot(lag_sec, task_mean_r(:, tk), '-', ...
        'Color', [C_LINE(tk, :) A_NS], 'LineWidth', LW_NS);

    segs = find_segments(task_sig0(:, tk));
    for s = 1:size(segs, 1)
        i1 = segs(s,1);   i2 = segs(s,2);
        plot(lag_sec(i1:i2), task_mean_r(i1:i2, tk), '-', ...
            'Color', C_LINE(tk, :), 'LineWidth', LW_SIG, 'HandleVisibility', 'off');
    end
end

for tk = 1:n_task
    plot(task_peak_lag(tk), task_peak_r(tk), 'o', 'MarkerSize', 6.5, ...
        'MarkerFaceColor', C_PEAK, 'MarkerEdgeColor', C_PEAK, ...
        'HandleVisibility', 'off');
end

h_sig = plot(nan, nan, '-', 'Color', [0.35 0.35 0.35], 'LineWidth', LW_SIG);
h_pk  = plot(nan, nan, 'o', 'MarkerSize', 6.5, 'LineStyle', 'none', ...
    'MarkerFaceColor', C_PEAK, 'MarkerEdgeColor', C_PEAK);

leg_str = [arrayfun(@(tk) sprintf('%s (N=%d)', TASK_NAMES{tk}, task_n(tk)), ...
           1:n_task, 'UniformOutput', false), ...
           {sprintf('r > 0 (p < %.2f)', ALPHA), 'Peak'}];
legend([h, h_sig, h_pk], leg_str, ...
    'Location', 'northwest', 'FontSize', 10, 'Box', 'on');

xlim([-MAX_LAG_SEC MAX_LAG_SEC]);
ya = [ci_lo_r(:); ci_hi_r(:)];
ylim([min(ya, [], 'omitnan') - 0.04, max(ya, [], 'omitnan') + 0.06]);

xlabel('Lag (s)   \leftarrow CMC leads | SMR leads \rightarrow', 'FontSize', 12);
ylabel('Cross-correlation (r)', 'FontSize', 12);
grid on; box on;
set(gca, 'FontSize', 11, 'LineWidth', 1.0, 'TickDir', 'out', 'Layer', 'top');

print(fig1, fullfile(OUTPUT_DIR, [OUT_PREFIX '_lagwise_curves']), '-dpng', '-r300');
savefig(fig1, fullfile(OUTPUT_DIR, [OUT_PREFIX '_lagwise_curves.fig']));

%% Figure 2: individual r at task-specific peak lags
BOX_W    = 0.55;
CAP_W    = 0.16;
JIT_W    = 0.26;
DOT_SIZE = 22;
LW       = 1.3;

fig2 = figure('Color', 'w', 'Position', [140 140 460 560]);
ax2  = axes('Parent', fig2); hold(ax2, 'on');
rng(PERM_SEED);

for tk = 1:n_task
    y = Rw(:, tk);
    y = y(~isnan(y));

    q1 = prctile(y, 25);   q2 = median(y);   q3 = prctile(y, 75);
    iqr_ = q3 - q1;
    lo_w = min(y(y >= q1 - 1.5*iqr_));
    hi_w = max(y(y <= q3 + 1.5*iqr_));

    xl = tk - BOX_W/2;   xr = tk + BOX_W/2;

    plot([tk tk], [hi_w q3], 'k-', 'LineWidth', LW);
    plot([tk tk], [q1 lo_w], 'k-', 'LineWidth', LW);
    plot([tk-CAP_W/2 tk+CAP_W/2], [hi_w hi_w], 'k-', 'LineWidth', LW);
    plot([tk-CAP_W/2 tk+CAP_W/2], [lo_w lo_w], 'k-', 'LineWidth', LW);

    patch([xl xr xr xl], [q1 q1 q3 q3], C_TASK(tk, :), ...
        'FaceAlpha', 1, 'EdgeColor', 'k', 'LineWidth', LW);
    plot([xl xr], [q2 q2], 'k-', 'LineWidth', LW);

    jit = (rand(numel(y), 1) - 0.5) * JIT_W;
    scatter(tk + jit, y, DOT_SIZE, 'k', 'filled', 'MarkerFaceAlpha', 0.9);
end

y_max  = max(Rw(:));
y_min  = min(Rw(:));
y_span = max(y_max - y_min, 0.1);
y_step = 0.11 * y_span;

lvl = 0;
for pp = 1:n_pair
    if p_dunn(pp) >= ALPHA, continue; end
    lvl = lvl + 1;
    xa = pair_list(pp,1);   xb = pair_list(pp,2);
    yb = y_max + y_step * lvl;
    plot([xa xb], [yb yb], 'k-', 'LineWidth', 1.1);
    text(mean([xa xb]), yb + 0.02*y_span, strtrim(star_str(p_dunn(pp))), ...
        'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
end

xlim([0.4, n_task + 0.6]);
xticks(1:n_task);
xticklabels(TASK_NAMES);
ylim([min(0, y_min - 0.06*y_span), y_max + y_step*(lvl + 1.1)]);
ylabel('Cross-correlation r at peak lag', 'FontSize', 13);

set(ax2, 'FontSize', 13, 'FontName', 'Arial', 'LineWidth', LW, ...
    'TickDir', 'out', 'TickLength', [0.018 0.018], ...
    'Box', 'off', 'XGrid', 'off', 'YGrid', 'off', 'Layer', 'top');

print(fig2, fullfile(OUTPUT_DIR, [OUT_PREFIX '_peaklag_comparison']), '-dpng', '-r300');
savefig(fig2, fullfile(OUTPUT_DIR, [OUT_PREFIX '_peaklag_comparison.fig']));

%% Figure 3: individual peak-lag distribution
fig3 = figure('Color', 'w', 'Position', [220 140 460 480]);
ax3  = axes('Parent', fig3); hold(ax3, 'on');
rng(PERM_SEED);

for tk = 1:n_task
    y = Lw(:, tk);
    q1 = prctile(y, 25);   q2 = median(y);   q3 = prctile(y, 75);
    iqr_ = q3 - q1;
    lo_w = min(y(y >= q1 - 1.5*iqr_));
    hi_w = max(y(y <= q3 + 1.5*iqr_));

    xl = tk - BOX_W/2;   xr = tk + BOX_W/2;

    plot([tk tk], [hi_w q3], 'k-', 'LineWidth', LW);
    plot([tk tk], [q1 lo_w], 'k-', 'LineWidth', LW);
    plot([tk-CAP_W/2 tk+CAP_W/2], [hi_w hi_w], 'k-', 'LineWidth', LW);
    plot([tk-CAP_W/2 tk+CAP_W/2], [lo_w lo_w], 'k-', 'LineWidth', LW);
    patch([xl xr xr xl], [q1 q1 q3 q3], C_TASK(tk, :), ...
        'FaceAlpha', 1, 'EdgeColor', 'k', 'LineWidth', LW);
    plot([xl xr], [q2 q2], 'k-', 'LineWidth', LW);

    jit = (rand(numel(y), 1) - 0.5) * JIT_W;
    scatter(tk + jit, y, DOT_SIZE, 'k', 'filled', 'MarkerFaceAlpha', 0.9);

    plot([xl xr], [task_peak_lag(tk) task_peak_lag(tk)], '-', ...
        'Color', C_PEAK, 'LineWidth', 1.8);
end

yline(0, 'k--', 'LineWidth', 0.8, 'Alpha', 0.4);

xlim([0.4, n_task + 0.6]);
xticks(1:n_task);
xticklabels(TASK_NAMES);
ylim([-MAX_LAG_SEC MAX_LAG_SEC]);
ylabel('Individual peak lag (s)', 'FontSize', 13);

set(ax3, 'FontSize', 13, 'FontName', 'Arial', 'LineWidth', LW, ...
    'TickDir', 'out', 'TickLength', [0.018 0.018], ...
    'Box', 'off', 'XGrid', 'off', 'YGrid', 'off', 'Layer', 'top');

print(fig3, fullfile(OUTPUT_DIR, [OUT_PREFIX '_individual_peak_lag']), '-dpng', '-r300');
savefig(fig3, fullfile(OUTPUT_DIR, [OUT_PREFIX '_individual_peak_lag.fig']));

fprintf('\nResults saved to:\n%s\n', OUTPUT_DIR);
fprintf('  %s_lagwise.csv\n',                    OUT_PREFIX);
fprintf('  %s_cluster_permutation.csv\n',        OUT_PREFIX);
fprintf('  %s_dunn_at_task_peak_lag.csv\n',      OUT_PREFIX);
fprintf('  %s_observation_level.csv\n',          OUT_PREFIX);
fprintf('  %s_results.mat\n',                    OUT_PREFIX);
fprintf('  %s_lagwise_curves.png/.fig\n',        OUT_PREFIX);
fprintf('  %s_peaklag_comparison.png/.fig\n',    OUT_PREFIX);
fprintf('  %s_individual_peak_lag.png/.fig\n',   OUT_PREFIX);

%% Local functions
function S = friedman_dunn(X, names)
% Friedman test with Dunn's post-hoc multiple comparisons.
%   X     : [n x k], rows = participants, columns = conditions, no missing data.
%   names : 1 x k condition names.
%
% Friedman: within-participant ranking with tie correction.
% Dunn's: two-tailed pairwise test with Bonferroni adjustment.

    [n, k] = size(X);

    % Rank conditions within each participant; tiedrank handles ties.
    R = zeros(n, k);
    for i = 1:n
        R(i, :) = tiedrank(X(i, :));
    end
    rank_sum = sum(R, 1);

    % Friedman statistic
    chi2 = 12 / (n * k * (k + 1)) * sum(rank_sum.^2) - 3 * n * (k + 1);

    % Tie correction
    tie_term = 0;
    for i = 1:n
        [~, ~, ic] = unique(X(i, :));
        cnt = accumarray(ic, 1);
        tie_term = tie_term + sum(cnt.^3 - cnt);
    end
    denom = 1 - tie_term / (n * k * (k^2 - 1));
    if denom > 0, chi2 = chi2 / denom; end

    df = k - 1;
    p  = 1 - chi2cdf(chi2, df);

    % Kendall's W
    W = chi2 / (n * (k - 1));

    % Dunn's post-hoc comparisons
    pair_list = nchoosek(1:k, 2);
    n_pair    = size(pair_list, 1);

    SE = sqrt(n * k * (k + 1) / 6);

    pair_A    = cell(n_pair, 1);
    pair_B    = cell(n_pair, 1);
    rank_diff = nan(n_pair, 1);
    z         = nan(n_pair, 1);
    p_unadj   = nan(n_pair, 1);
    p_adj     = nan(n_pair, 1);

    for pp = 1:n_pair
        a = pair_list(pp, 1);
        b = pair_list(pp, 2);

        pair_A{pp}    = names{a};
        pair_B{pp}    = names{b};
        rank_diff(pp) = rank_sum(a) - rank_sum(b);
        z(pp)         = abs(rank_diff(pp)) / SE;
        p_unadj(pp)   = 2 * (1 - normcdf(z(pp)));
        p_adj(pp)     = min(p_unadj(pp) * n_pair, 1);
    end

    T_pair = table(pair_A, pair_B, rank_diff, z, p_unadj, p_adj, ...
        p_adj < 0.05, ...
        'VariableNames', {'TaskA', 'TaskB', 'RankSumDiff', 'z', ...
                          'p_unadjusted', 'p_adjusted', 'Significant'});

    S = struct('n', n, 'k', k, 'df', df, 'chi2', chi2, 'p', p, 'W', W, ...
               'rank_sum', rank_sum, 'rank_matrix', R, ...
               'pair_list', pair_list, 'n_pair', n_pair, ...
               'pair_A', {pair_A}, 'pair_B', {pair_B}, ...
               'rank_diff', rank_diff, 'z', z, ...
               'p_unadj', p_unadj, 'p_adj', p_adj, 'T_pair', T_pair);
end

function [masses, segs] = cluster_mass_right(t_vec, tcrit)
% One-tailed positive supra-threshold clusters and cluster mass (sum of t).
    mask = t_vec > tcrit;
    idx  = find(mask);
    if isempty(idx), masses = []; segs = []; return; end

    bp    = find(diff(idx) > 1);
    seg_s = [idx(1); idx(bp + 1)];
    seg_e = [idx(bp); idx(end)];
    segs  = [seg_s, seg_e];

    masses = arrayfun(@(a, b) sum(t_vec(a:b)), seg_s, seg_e);
end

function segs = find_segments(mask)
    mask = logical(mask(:));
    idx = find(mask);
    if isempty(idx), segs = []; return; end
    bp = find(diff(idx) > 1);
    seg_s = [idx(1); idx(bp + 1)];
    seg_e = [idx(bp); idx(end)];
    segs = [seg_s, seg_e];
end

function s = star_str(p)
    if     p < 0.001, s = ' ***';
    elseif p < 0.01,  s = ' **';
    elseif p < 0.05,  s = ' *';
    else,             s = '';
    end
end