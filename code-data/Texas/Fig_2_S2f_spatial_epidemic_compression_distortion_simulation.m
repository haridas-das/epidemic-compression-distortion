%% =========================================================================
%  Fig_2_S2f_epidemic_compression_distortion_simulation.m
%
%  County-network SIR simulations for Fig. 2 (main text) and the network-
%  density sensitivity analysis (Supplementary Materials) of:
%    "Compression versus distortion: choosing the spatial scale of
%     aggregation in epidemic models"
%
%  For one state (Texas, Colorado, or Massachusetts) the script:
%    1. loads the daily 2020 county-to-county flux matrices,
%    2. thresholds every day at phi^th = 0.02, resets the diagonal so rows
%       sum to one, and selects the least connected day as a static network,
%    3. simulates the metapopulation SIR-network model on that network,
%    4. aggregates prevalence to the county, MSA, and state levels,
%    5. repeats the simulation retaining the strongest 5/10/15/100% of links
%       of the same day, and for a range of fixed thresholds, reporting the
%       mean absolute optimal lag (Eq. 7.2) of each aggregate.
%
%  HOW TO RUN
%    1. Set STATE in Section 0 ('TX', 'CO', or 'MA').
%    2. Run the script from that state's data folder, which must contain
%       county_name_list.m, pop_vector_origin_only.csv, and the
%       <County>_Mobility_Flux.csv files.
%    3. Outputs (figures and CSV tables) are written to ./result_figures.
%
%  REQUIRES
%    Functions: county_name_list.m, TDFluxMatrix.m,
%               InitialValueMobilitySIRNetwork.m, SirGenNetworkModel.m
%    MATLAB R2018b or later (xline, islocalmax); R2020a or later for the
%    vector PDF export (exportgraphics). No toolboxes are required.
%    Tested with MATLAB R2024b.
% =========================================================================
clear; close all; clc;

%% 0. SETTINGS ============================================================
% >>>>>>>>>>>>>>>>>>>>>>>  CHOOSE THE STATE HERE  <<<<<<<<<<<<<<<<<<<<<<<<<<
STATE = 'TX';            % 'TX' = Texas, 'CO' = Colorado, 'MA' = Massachusetts

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% Only the county count, county names, first data day, and the top-five
% nodes are taken from county_name_list (of the state folder you run in).
% Its first output is the state FIPS code (e.g. 8), so the name is set below.
[~, m, county_names, time_start, ~, top5_nodes] = county_name_list();

%% 0b. STATE-SPECIFIC SETTINGS -------------------------------------------
%  n_counties  : expected number of counties (used to check the data folder)
%  msa_counties: counties summed for the MSA curve (solid red line); in each
%                state this is the largest metro area (TX: Dallas-Fort Worth,
%                CO: Denver, MA: Boston), as in the original Fig. 2 scripts
%  stay_home   : statewide stay-at-home period, ONLY shaded in the mobility
%                figure; it does not affect the analysis. [] = no shading.
%  flux_pair   : {origin, destination} for the county-pair flux figure
switch STATE
    case 'CO'
        state_name   = 'Colorado';
        n_counties   = 64;
        % Denver metro area, 7 counties -- the list used in the original
        % Colorado Fig. 2 script (..._MSA_COfi.m)
        msa_counties = {'Adams County','Arapahoe County','Boulder County', ...
                        'Broomfield County','Denver County','Douglas County', ...
                        'Jefferson County'};
        stay_home    = [datetime(2020,3,26) datetime(2020,4,26)];   % statewide stay-at-home order
        flux_pair    = {'Adams County','Denver County'};
    case 'TX'
        state_name   = 'Texas';
        n_counties   = 254;
        % Dallas-Fort Worth-Arlington MSA, 13 counties -- the list used in the
        % original Texas Fig. 2 script (..._MSA_TXfi.m). NOTE: this is the
        % state's largest metro area, NOT the metro area of the seed county
        % (Harris County is in the Houston MSA).
        msa_counties = {'Collin County','Dallas County','Denton County','Ellis County', ...
                        'Hunt County','Kaufman County','Rockwall County','Hood County', ...
                        'Johnson County','Parker County','Somervell County', ...
                        'Tarrant County','Wise County'};
        stay_home    = [];      % no conventional statewide stay-at-home order
        flux_pair    = {'Fort Bend County','Harris County'};
    case 'MA'
        state_name   = 'Massachusetts';
        n_counties   = 14;
        % Boston metro area, 4 counties -- the list used in the original
        % Massachusetts Fig. 2 panel (Plymouth County NOT included; the
        % official OMB Boston-Cambridge-Newton MSA also contains Plymouth).
        % Order does not matter: the MSA curve is the sum over these counties.
        msa_counties = {'Middlesex County', 'Essex County', 'Suffolk County', 'Norfolk County'};
        stay_home    = [];      % no shading
        flux_pair    = {'Middlesex County','Suffolk County'};
    otherwise
        error('STATE must be ''TX'', ''CO'' or ''MA''.');
end
state_abbr = STATE;                       % used in output file names

% Safety check: are we in the right state's data folder?
if m ~= n_counties
    error(['STATE = ''%s'' expects %d counties, but county_name_list returned %d.\n' ...
           'Run the script from the %s data folder, or change STATE.'], ...
           STATE, n_counties, m, state_name);
end

% --- Epidemic parameters -------------------------------------------------
beta0 = 0.20;            % transmission rate (day^-1)
gamma = 0.10;            % recovery rate (day^-1)  -> local R0 = 2
zeta  = 1;               % multiplier for the seed county's beta (1 = homogeneous)
Nv    = 1;               % number of initially infectious individuals
seed  = top5_nodes(1);   % seed county: most populous (Harris / Denver / Middlesex)

% --- Mobility network ----------------------------------------------------
mob_day_end = 366;       % last day of 2020 flux data loaded (2020 is a leap year)
phi_th      = 0.02;      % sparsification threshold phi^th
zero_tol    = 1e-12;     % entries below this count as zero when measuring sparsity

% Diagonal after thresholding:
%   true  -> (DEFAULT) phi_ii = 1 - sum_{j~=i} phi_ij, so every row of the
%            simulated network sums to one (row-stochastic, as in Methods 7.1);
%            mobility removed by the threshold counts as staying home.
%   false -> reproduces the ORIGINAL Fig. 2 runs: the thresholded matrix is
%            used as is, so rows sum to LESS than one. (The original scripts
%            computed phi_ii = 1 - sum_{j~=i} phi_ij, but the next line,
%            "matrix_phi = StaticTDNetwork", overwrote it.)
reset_diagonal = true;

% --- Calendar alignment of the flux files --------------------------------
% Row 1 of every <County>_Mobility_Flux.csv is taken to be 1 Jan 2020.
% Section 2b checks this with the data (weekly cycle, March 2020 collapse).
data_date0 = datetime(2020,1,1);
% (The stay-at-home period used for shading is set per state in Section 0b.)

% --- Simulation horizon --------------------------------------------------
dt      = 0.1;           % output step (days); ode45 chooses its own internal steps
T_end   = 500;           % simulated days (network is static, so no calendar limit)
tspan   = 1:dt:T_end;
x_panel = [50 450];      % time window shown in the Fig. 2-style panels

% --- Density sweep (Section 5) -------------------------------------------
run_density_sweep = true;
sweep_fixed = [0.05 0.10 0.15 1.00];      % shares of nonzero links kept (panels besides phi^th)
                                          % the phi^th panel is inserted at its true
                                          % density automatically (Section 5)
lag_max = 90;                             % lag search bound (days), as in Eq. 7.2

% --- Solver tolerances (used for EVERY simulation in this script) --------
% true : relative 1e-7, absolute 1e-9, as in the canonical simulations (Methods 7.1)
% false: MATLAB ode45 defaults (relative 1e-3, absolute 1e-6), as in the original runs
tight_tol = true;
if tight_tol
    ode_opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);
else
    ode_opts = odeset();
end

% --- Extra analyses -------------------------------------------------------
run_lcc_check       = true;                                  % Section 2d
run_threshold_sweep = true;                                  % Section 6
th_list = [0 0.002 0.005 0.01 0.015 0.02 0.03 0.05 0.08];    % thresholds for Section 6 (0 = all links)

% --- Output folders and file tag -----------------------------------------
out_dir = fullfile(pwd, 'result_figures');
if ~exist(fullfile(out_dir,'fig'), 'dir'), mkdir(fullfile(out_dir,'fig')); end
if ~exist(fullfile(out_dir,'png'), 'dir'), mkdir(fullfile(out_dir,'png')); end
net_tag = ['StaticTDNetwork_With_threshold_' num2str(zero_tol) ...
           '_inter_mob_threshold_' num2str(phi_th)];
if reset_diagonal                         % keep the two versions' outputs apart
    net_tag    = [net_tag '_rowstoch'];
    state_abbr = [state_abbr '_rowstoch'];
end
fig_tag = sprintf('%s_seed%d_beta%g_gamma%g', net_tag, seed, beta0, gamma);


%% 1. POPULATION AND INITIAL CONDITIONS ==================================
pop   = csvread('pop_vector_origin_only.csv', 1, 0);   % skip header row
vec_N = pop(:, 2);                                     % county populations

[vec_S0, vec_I0, vec_N, y0] = InitialValueMobilitySIRNetwork(m, vec_N, seed, Nv);

vec_beta       = beta0 * ones(m, 1);   % same beta in every county ...
vec_beta(seed) = zeta * beta0;         % ... unless zeta ~= 1 for the seed


%% 2. STATIC MOBILITY NETWORK ============================================
% (a) Load the daily 2020 flux matrices: Phi_raw(i,j,k) = fraction of
%     county i's residents who travel to county j on data day k.
Phi_raw = TDFluxMatrix(m, m, 1, time_start, mob_day_end, county_names);
n_days  = size(Phi_raw, 3);

% Check: rows of the RAW daily matrices should sum to one. The diagonal
% (staying home) was already set to 1 - (off-diagonal sum) for all 366 days
% when the <County>_Mobility_Flux.csv files were prepared.
rs_raw = squeeze(sum(Phi_raw, 2));                 % m x n_days
fprintf('Raw daily matrices, row sums: min %.4f, max %.4f\n', min(rs_raw(:)), max(rs_raw(:)));

% (b) Sparsify every day: entries below phi^th are set to zero.
Phi_thr = Phi_raw;
Phi_thr(Phi_thr < phi_th) = 0;
rs_cut = squeeze(sum(Phi_thr, 2));                 % row sums right after thresholding
fprintf('Thresholded daily matrices BEFORE diagonal reset, row sums: min %.4f, max %.4f\n', ...
        min(rs_cut(:)), max(rs_cut(:)));

% (b2) Adjust the diagonal on EVERY day after thresholding (reset_diagonal = true):
%      phi_ii(t) = 1 - sum_{j~=i} phi_ij(t), so each daily matrix stays
%      row-stochastic. Mobility removed by the threshold counts as staying home.
fprintf('Smallest raw diagonal entry (any county, any day): %.4f\n', ...
        min(arrayfun(@(k) min(diag(Phi_raw(:,:,k))), 1:n_days)));
if reset_diagonal
    for k = 1:n_days
        Phi_thr(:,:,k) = make_row_stochastic(Phi_thr(:,:,k));
    end
    rs_thr = squeeze(sum(Phi_thr, 2));             % m x n_days, should all be 1
    fprintf('Thresholded daily matrices after diagonal reset, row sums: min %.4f, max %.4f (all %d days)\n', ...
            min(rs_thr(:)), max(rs_thr(:)), n_days);
end

% (c) Pick the sparsest thresholded day (largest fraction of zero entries).
zero_frac = zeros(1, n_days);
for k = 1:n_days
    zero_frac(k) = nnz(Phi_thr(:,:,k) < zero_tol) / (m*m);
end
[~, k_sparse] = max(zero_frac);

% (d) This single matrix is the static network used for the whole run.
%     With reset_diagonal = true (default) its diagonal was already adjusted
%     in step (b2), so rows sum to one; with false its rows sum to less than
%     one (original Fig. 2 runs).
Phi_static = Phi_thr(:,:,k_sparse);
rs = sum(Phi_static, 2);
fprintf('Row sums of the simulated network: min %.4f, mean %.4f, max %.4f\n', ...
        min(rs), mean(rs), max(rs));

% --- Share of links kept by phi^th on the selected day ----------------
off     = ~eye(m);                                 % off-diagonal mask
Phi_day = Phi_raw(:,:,k_sparse);                   % selected day, before thresholding
n_links = nnz(Phi_day(off) > 0);                   % nonzero inter-county links
n_kept  = nnz(Phi_static(off) > 0);                % links surviving phi^th
fprintf('Links kept at phi_th = %.2f: %.1f%% of nonzero links, %.1f%% of all m(m-1) pairs\n', ...
        phi_th, 100*n_kept/n_links, 100*n_kept/(m*(m-1)));


%% 2b. SNAPSHOT DATE ===============================================
% Data day k is row (time_start + k - 1) of the flux files, i.e. the date
% data_date0 + (time_start + k - 2) days.
dates     = data_date0 + days((time_start:time_start+n_days-1) - 1);
snap_date = dates(k_sparse);
n_tie     = nnz(abs(zero_frac - zero_frac(k_sparse)) < 1e-12);

fprintf('\n%s: selected snapshot = %s (%s), data day %d, %.2f%% zero entries\n', ...
        state_name, datestr(snap_date, 'dd mmm yyyy'), char(day(snap_date, 'name')), ...
        k_sparse, 100*zero_frac(k_sparse));
if n_tie > 1
    fprintf('  Note: %d days share this sparsity; MATLAB kept the earliest.\n', n_tie);
end

% --- The 10 sparsest days (is the choice clear-cut?) ---------------------
[z_sorted, order] = sort(zero_frac, 'descend');
n_show = min(10, n_days);
T_top  = table(dates(order(1:n_show))', string(day(dates(order(1:n_show)), 'name'))', ...
               100*z_sorted(1:n_show)', 'VariableNames', {'date','weekday','pct_zero_entries'});
disp(T_top)

% --- Daily mobility index -----------------------------------------------
% Average share of a county's residents traveling to another county that
% day (raw flux, before thresholding).
outflow  = squeeze(sum(Phi_raw .* off, 2));        % m x n_days
mob_idx  = mean(outflow, 1)';                      % n_days x 1
% Number of inter-county links (i ~= j) with phi_ij >= phi^th on each day.
% The diagonal is always kept, so the day with the most zero entries is the
% day with the FEWEST such links -- the selected snapshot.
n_links_day = squeeze(sum(sum((Phi_thr > 0) & off, 1), 2));   % n_days x 1

% --- Check 1: weekly cycle (lowest mobility should fall on Sunday) ------
wd     = weekday(dates)';                          % 1 = Sunday ... 7 = Saturday
wd_avg = arrayfun(@(k) mean(mob_idx(wd == k)), 1:7);
[~, wd_low] = min(wd_avg);
wd_names = {'Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'};
fprintf('Check 1 - lowest-mobility weekday: %s (expected: Sunday)\n', wd_names{wd_low});

% --- Check 2: mobility collapse in March 2020 ----------------------------
pre   = mean(mob_idx(dates >= datetime(2020,3,1) & dates <= datetime(2020,3,10)));
lock  = mean(mob_idx(dates >= datetime(2020,4,1) & dates <= datetime(2020,4,15)));
[~, kd] = min(diff(movmean(mob_idx, 7)));          % sharpest fall of the 7-day average
fprintf('Check 2 - mobility 1-10 Mar vs 1-15 Apr: %.3f vs %.3f (%.0f%% drop); sharpest fall: %s (expected: mid/late March)\n', ...
        pre, lock, 100*(1 - lock/pre), datestr(dates(kd+1), 'dd mmm yyyy'));

% --- Save the snapshot information --------------------------------------
T_snap = table(string(state_name), snap_date, string(day(snap_date,'name')), k_sparse, ...
               time_start + k_sparse - 1, 100*zero_frac(k_sparse), n_tie, ...
               100*n_kept/n_links, 100*n_kept/(m*(m-1)), 'VariableNames', ...
               {'state','snapshot_date','weekday','data_day','file_row','pct_zero_entries', ...
                'n_tied_days','links_pct_of_nonzero','links_pct_of_all_pairs'});
writetable(T_snap, fullfile(out_dir, sprintf('snapshot_date_%s.csv', state_abbr)));
writetable(T_top,  fullfile(out_dir, sprintf('sparsest_10_days_%s.csv', state_abbr)));

% --- SI figure: daily mobility with the selected day marked --------------
% Dates are plotted as datenum so the shaded stay-at-home band works in any
% MATLAB version (patch does not accept datetime on older releases).
xd = datenum(dates);  x0 = datenum(snap_date);
has_band = ~isempty(stay_home);
if has_band, xs = datenum(stay_home); end
hD = figure('Color','w', 'Units','pixels', 'Position',[100 100 1400 700]);
for pnl = 1:2
    subplot(2, 1, pnl); hold on
    if pnl == 1, yv = mob_idx;  ylab = {'Mean share of residents','leaving their county'};
    else,        yv = n_links_day; ylab = {'Inter-county links', sprintf('with \\phi_{ij} \\geq %.2f', phi_th)};
    end
    yl = [min(yv) max(yv)] + [-0.05 0.05]*(max(yv) - min(yv));
    if has_band
        patch([xs fliplr(xs)], yl([1 1 2 2]), [0.85 0.85 0.85], 'EdgeColor','none', 'FaceAlpha', 0.6);
    end
    plot(xd, yv, 'k-', 'LineWidth', 1.2);
    plot(xd, movmean(yv, 7), 'b-', 'LineWidth', 2.5);
    plot(x0, yv(k_sparse), 'ro', 'MarkerSize', 12, 'LineWidth', 2.5);
    xlim([xd(1) xd(end)]); ylim(yl); ylabel(ylab);
    datetick('x', 'mmm', 'keeplimits');
    style_axes(gca, 18); hold off
    if pnl == 1
        leg = {'Daily','7-day average', sprintf('Selected snapshot (%s)', datestr(snap_date,'dd mmm yyyy'))};
        if has_band, leg = [{'Stay-at-home period'}, leg]; end
        legend(leg, 'Location', 'northeast');
        title(sprintf('%s: daily inter-county mobility, 2020', state_name), 'FontWeight', 'normal');
    end
end
save_fig(hD, out_dir, sprintf('daily_mobility_snapshot_%s', state_abbr));


%% 2c. DAILY FLUX BETWEEN TWO COUNTIES (SI figure) ==================
% Time series of phi_ij(t) for the pair chosen in Section 0b, both directions:
% the fraction of residents of the origin county traveling to the destination.
io = county_index(county_names, flux_pair(1));
jd = county_index(county_names, flux_pair(2));
if isempty(io) || isempty(jd)
    warning('Flux pair not found in county list; pair figure skipped.');
else
    nm     = erase(string(flux_pair), " County");            % short names for labels
    series = {squeeze(Phi_raw(io, jd, :)), squeeze(Phi_raw(jd, io, :))};
    labels = {sprintf('%s \\rightarrow %s', nm(1), nm(2)), ...
              sprintf('%s \\rightarrow %s', nm(2), nm(1))};
    hQ = figure('Color','w', 'Units','pixels', 'Position',[100 100 1400 700]);
    for p = 1:2
        subplot(2, 1, p); hold on
        y = series{p};
        plot(dates, y, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);                 % daily value
        plot(dates, movmean(y, 7), 'b-', 'LineWidth', 2.5);                      % 7-day average
        plot(dates([1 end]), [phi_th phi_th], 'k--', 'LineWidth', 1.5);           % threshold
        plot(snap_date, y(k_sparse), 'ro', 'MarkerSize', 12, 'LineWidth', 2.5);  % selected day
        ylabel({'Fraction of origin','residents traveling'});
        title(labels{p}, 'FontWeight', 'normal');
        style_axes(gca, 18); hold off
        if p == 1
            legend({'Daily','7-day average', sprintf('\\phi^{th} = %.2f', phi_th), ...
                    sprintf('Selected day (%s)', datestr(snap_date,'dd mmm yyyy'))}, ...
                   'Location', 'northeast');
        end
    end
    save_fig(hQ, out_dir, sprintf('flux_%s_to_%s_%s', strrep(nm(1),' ','_'), ...
             strrep(nm(2),' ','_'), state_abbr));   % e.g. flux_Fort_Bend_to_Harris_TX
end


%% 2d. CONNECTIVITY VS. THRESHOLD (SI figure) =======================
% For the selected day: share of counties in the largest connected component
% (direction ignored) and share of nonzero inter-county links kept, as the
% threshold varies. A flat stretch of the first curve around phi^th supports
% the "connectivity plateau" statement in the Methods.
if run_lcc_check
    th_grid = logspace(-3, log10(0.2), 60);
    lcc = zeros(size(th_grid));  links_kept_grid = lcc;
    for q = 1:numel(th_grid)
        A = (Phi_day >= th_grid(q)) & off;                  % links kept at this threshold
        G = graph(A | A', 'omitselfloops');
        comp = conncomp(G);
        lcc(q) = max(histcounts(comp, 1:max(comp)+1)) / m;  % share of counties in largest component
        links_kept_grid(q) = 100 * nnz(A) / n_links;
    end
    hL = figure('Color','w');
    yyaxis left;  semilogx(th_grid, lcc, 'o-', 'LineWidth', 1.5);
    ylabel('Share of counties in largest component');
    yyaxis right; semilogx(th_grid, links_kept_grid, 's-', 'LineWidth', 1.5);
    ylabel('% of nonzero links kept');
    xline(phi_th, 'k--', sprintf('\\phi^{th} = %.2f', phi_th));
    xlabel('Threshold \phi^{th}');
    title(sprintf('%s, %s', state_name, datestr(snap_date, 'dd mmm yyyy')), 'FontWeight', 'normal');
    style_axes(gca, 16);
    save_fig(hL, out_dir, sprintf('lcc_vs_threshold_%s', state_abbr));
    writetable(table(th_grid(:), lcc(:), links_kept_grid(:), 'VariableNames', ...
               {'phi_th','share_in_largest_component','links_pct_of_nonzero'}), ...
               fullfile(out_dir, sprintf('lcc_vs_threshold_%s.csv', state_abbr)));
end


%% 3. SIMULATION ==========================================================
[t, Y] = ode45(@(t,y) SirGenNetworkModel(m, Phi_static, y, vec_beta, ...
               vec_S0, vec_I0, vec_N, gamma, zeta), tspan, y0, ode_opts);

I       = Y(:, m+1:2*m);                           % prevalence, one column per county
msa_idx = county_index(county_names, msa_counties);
I_state = sum(I, 2);                               % state aggregate
I_msa   = sum(I(:, msa_idx), 2);                   % MSA aggregate

% Fig. 2 panel -- counties (blue), MSA (solid red), state (dashed red)
h = figure('Color','w'); hold on
Z = minmax_cols(I);
others = setdiff(1:m, top5_nodes);
plot(t, Z(:, others), 'Color', [0 0 1 0.2], 'LineWidth', 0.8);
for k = 1:numel(top5_nodes)
    plot(t, Z(:, top5_nodes(k)), 'Color', [0 0 1], 'LineWidth', 5/k);
end
plot(t, minmax_cols(I_state), 'r--', 'LineWidth', 7);
plot(t, minmax_cols(I_msa),   'r-',  'LineWidth', 5);
xlim(x_panel); xticks(0:50:450);
style_axes(gca, 30); hold off
save_fig(h, out_dir, ['Fig_MSA_state_county_minmax_' fig_tag]);


%% 4. CSV OUTPUTS ========================================================

% Full prevalence time series: Time, Node_1, ..., Node_m
T_dyn = array2table([t, I], 'VariableNames', ...
        [{'Time'}, arrayfun(@(k) sprintf('Node_%d', k), 1:m, 'UniformOutput', false)]);
writetable(T_dyn, fullfile(out_dir, sprintf( ...
    'disease_dynamics_for_%s_I0_%g_gamma_%g_beta0_%g.csv', net_tag, seed, gamma, beta0)));

% Transmission rate per county
T_beta = table((1:m)', vec_beta, 'VariableNames', {'County_ID','Beta'});
writetable(T_beta, fullfile(out_dir, sprintf( ...
    'vec_beta_for_%s_I0_%g_gamma_%g_beta0_%g.csv', net_tag, seed, gamma, beta0)));


%% 5. DENSITY SWEEP (SI figure) ======================================
% Same model, seed, parameters and snapshot day as Fig. 2. Only the number
% of retained inter-county links changes:
%   strongest 5%  |  phi^th = 0.02 (Fig. 2)  |  10%  |  15%  |  all links
if run_density_sweep
    w       = Phi_day(off & Phi_day > 0);     % weights of all nonzero links that day
    day_idx = 1:round(1/dt):numel(t);         % daily samples (t = 1, 2, ..., T_end)

    % Place the phi^th (Fig. 2) panel at its true density among the fixed
    % shares, so panels always run from sparsest to densest in every state.
    share_th   = n_kept / n_links;                         % e.g. 0.071 for Colorado
    pos        = find(sweep_fixed > share_th, 1);          % first fixed share above it
    sweep_keep = [sweep_fixed(1:pos-1), NaN, sweep_fixed(pos:end)];
    fprintf('Density sweep order (NaN = phi^th panel, %.1f%% of links): %s\n', ...
            100*share_th, mat2str(sweep_keep));
    R  = table();
    hS = figure('Color','w', 'Units','pixels', 'Position',[50 50 2400 520]);

    for p = 1:numel(sweep_keep)
        % --- network for this panel -----------------------------------
        if isnan(sweep_keep(p))
            M  = Phi_static;                  % exactly the Fig. 2 network
            Ip = I;                           % reuse the Fig. 2 run
        else
            % Keep the strongest share of links by ranking the weights directly.
            ws     = sort(w, 'descend');                          % strongest first
            n_keep = max(1, round(sweep_keep(p) * numel(ws)));    % e.g. 5% of 1,199 = 60 links
            thr    = ws(n_keep);                                  % weight of weakest kept link
            M      = Phi_day;
            M(off & M < thr) = 0;
            fprintf('  keep %3.0f%%: cut-off phi = %.4g, links kept = %d of %d\n', ...
                    100*sweep_keep(p), thr, nnz(M(off) > 0), n_links);
            if reset_diagonal, M = make_row_stochastic(M); end   % same rule as Fig. 2
            [~, Yp] = ode45(@(t,y) SirGenNetworkModel(m, M, y, vec_beta, ...
                            vec_S0, vec_I0, vec_N, gamma, zeta), tspan, y0, ode_opts);
            Ip = Yp(:, m+1:2*m);
        end

        % Share of links kept (the definition used in the paper):
        % nonzero inter-county links of M / nonzero inter-county links of the
        % original matrix of the selected day. The diagonal is not a link.
        pct_links = 100 * nnz(M(off) > 0) / n_links;
        if isnan(sweep_keep(p)), cutoff = phi_th; else, cutoff = thr; end
        % ------------------------------------------------------------------

        % --- metrics (daily samples, Eq. 7.2) -------------------------
        Id  = Ip(day_idx, :);
        td  = t(day_idx);
        Is  = sum(Id, 2);
        Im  = sum(Id(:, msa_idx), 2);
        hit = find(max(Id, [], 1) > 1);       % counties with an actual outbreak
        [~, ip] = max(Id(:, hit), [], 1);
        lagS = mean_abs_lag(Is, Id, hit, lag_max);
        lagM = mean_abs_lag(Im, Id, intersect(hit, msa_idx), lag_max);
        % Lag including unreached counties: a county the epidemic never reaches
        % is given the maximum lag (lag_max), i.e. treated as maximally distorted.
        hitM     = intersect(hit, msa_idx);
        lagS_all = (lagS*numel(hit)  + lag_max*(m - numel(hit)))              / m;
        if isempty(hitM)                  % MSA never reached (e.g. Dallas in a sparse Texas network)
            lagM_all = lag_max;
        else
            lagM_all = (lagM*numel(hitM) + lag_max*(numel(msa_idx) - numel(hitM))) / numel(msa_idx);
        end
        npk  = nnz(islocalmax(minmax_cols(Is), 'MinProminence', 0.05));

        % --- panel title and table label ------------------------------
        if isnan(sweep_keep(p))
            ttl = sprintf('$\\phi^{th}=%.2f$ (%.1f\\%% of links, Fig. 2)', phi_th, pct_links);
            lab = sprintf("phi_th = %.2f", phi_th);
        elseif sweep_keep(p) == 1
            ttl = 'All links (no sparsification)';
            lab = "all links";
        else
            ttl = sprintf('Strongest %.0f\\%% of links', 100*sweep_keep(p));
            lab = sprintf("strongest %.0f%%", 100*sweep_keep(p));
        end

        R = [R; table(lab, cutoff, nnz(M(off) > 0), pct_links, 100*nnz(M(off) > 0)/(m*(m-1)), ...
                      numel(hit), std(td(ip)), lagS, lagM, lagS_all, lagM_all, npk, 'VariableNames', ...
                      {'network','cutoff_phi','links_kept','links_pct_of_nonzero', ...
                       'links_pct_of_all_pairs','counties_infected','sd_peak_days', ...
                       'state_lag_days','MSA_lag_days','state_lag_incl_unreached_days', ...
                       'MSA_lag_incl_unreached_days','state_peaks'})];

        % --- panel ----------------------------------------------------
        subplot(1, numel(sweep_keep), p); hold on
        plot(td, minmax_cols(Id(:, hit)), 'Color', [0 0 1 0.2], 'LineWidth', 0.8, ...
             'HandleVisibility', 'off');
        hC = plot(nan, nan, 'Color', [0.4 0.4 1], 'LineWidth', 1.5);   % legend entry
        hM = plot(td, minmax_cols(Im), 'r-',  'LineWidth', 4);
        hT = plot(td, minmax_cols(Is), 'r--', 'LineWidth', 4);
        xlim(x_panel); ylim([0 1]); xticks(100:100:400);
        style_axes(gca, 16);
        title(ttl, 'Interpreter', 'latex', 'FontSize', 18);
        text(0.97, 0.80, sprintf('State lag: %.1f d (%.1f d)\nMSA lag: %.1f d\nInfected: %d/%d', ...
             lagS, lagS_all, lagM, numel(hit), m), 'Units', 'normalized', ...
             'HorizontalAlignment', 'right', 'FontSize', 14, 'FontName', 'Times');
        xlabel('Time (days)');
        if p == 1, ylabel('Normalized cases'); end
        if p == numel(sweep_keep)
            legend([hC hM hT], {'County','MSA','State'}, 'Location', 'east');
        end
        hold off
    end

    disp(R)
    writetable(R, fullfile(out_dir, sprintf('density_sweep_%s.csv', state_abbr)));
    save_fig(hS, out_dir, sprintf('density_sweep_%s', state_abbr));
    try   % vector PDF for the SI (MATLAB R2020a or later)
        exportgraphics(hS, fullfile(out_dir, sprintf('density_sweep_%s.pdf', state_abbr)), ...
                       'ContentType', 'vector');
    catch
        warning('exportgraphics unavailable; use the .png/.fig versions.');
    end
end


%% 6. DISTORTION VS. THRESHOLD phi^th (SI figure) ==================
% Same day, seed, and parameters; links below each fixed threshold removed.
% Expectation: distortion of the state aggregate grows with phi^th. The lag
% including unreached counties should rise steadily; the lag over reached
% counties alone can fall again once the network breaks apart.
if run_threshold_sweep
    day_idx = 1:round(1/dt):numel(t);
    T_th = table();
    for q = 1:numel(th_list)
        M = Phi_day;  M(off & M < th_list(q)) = 0;
        if reset_diagonal, M = make_row_stochastic(M); end
        [~, Yq] = ode45(@(t,y) SirGenNetworkModel(m, M, y, vec_beta, ...
                        vec_S0, vec_I0, vec_N, gamma, zeta), tspan, y0, ode_opts);
        Iq  = Yq(day_idx, m+1:2*m);
        Isq = sum(Iq, 2);  Imq = sum(Iq(:, msa_idx), 2);
        hq  = find(max(Iq, [], 1) > 1);
        hqM = intersect(hq, msa_idx);
        lS  = mean_abs_lag(Isq, Iq, hq,  lag_max);
        lM  = mean_abs_lag(Imq, Iq, hqM, lag_max);
        T_th = [T_th; table(th_list(q), nnz(M(off) > 0), 100*nnz(M(off) > 0)/n_links, numel(hq), ...
                    lS, (lS*numel(hq) + lag_max*(m - numel(hq)))/m, lM, ...
                    'VariableNames', {'phi_th','links_kept','links_pct_of_nonzero', ...
                    'counties_reached','state_lag_days','state_lag_incl_unreached_days', ...
                    'MSA_lag_days'})]; %#ok<AGROW>
        fprintf('  phi_th = %.3f: %d links, %d/%d counties reached, state lag %.1f d (%.1f d incl. unreached)\n', ...
                th_list(q), nnz(M(off) > 0), numel(hq), m, lS, T_th.state_lag_incl_unreached_days(end));
    end
    disp(T_th)
    writetable(T_th, fullfile(out_dir, sprintf('distortion_vs_threshold_%s.csv', state_abbr)));

    hT = figure('Color','w'); hold on
    yyaxis left
    plot(T_th.phi_th, T_th.state_lag_incl_unreached_days, 'o-', 'LineWidth', 2);
    plot(T_th.phi_th, T_th.state_lag_days, 's--', 'LineWidth', 2);
    ylabel('State lag (days)');
    yyaxis right
    plot(T_th.phi_th, T_th.counties_reached, 'd:', 'LineWidth', 2);
    ylabel('Counties reached');
    xline(phi_th, 'k--', 'Fig. 2');
    xlabel('Threshold \phi^{th}');
    legend({'State lag, unreached = max lag', 'State lag, reached counties only', ...
            'Counties reached'}, 'Location', 'best');
    title(sprintf('%s: distortion vs. threshold', state_name), 'FontWeight', 'normal');
    style_axes(gca, 16); hold off
    save_fig(hT, out_dir, sprintf('distortion_vs_threshold_%s', state_abbr));
end


%% =========================================================================
%  LOCAL FUNCTIONS (must stay at the end of the file)
%  =========================================================================
function Z = minmax_cols(X)
% Min-Max normalize each column to [0,1]; constant columns become zeros.
lo = min(X, [], 1);  rng = max(X, [], 1) - lo;
Z  = zeros(size(X));
ok = rng > 0;
Z(:, ok) = (X(:, ok) - lo(ok)) ./ rng(ok);
end

function idx = county_index(names, wanted)
% Indices of the requested counties in the model's county list.
idx = [];
for k = 1:numel(wanted)
    j = find(strcmp(names, wanted{k}));
    if isempty(j)
        warning('County "%s" not found in county list.', wanted{k});
    else
        idx(end+1) = j; %#ok<AGROW>
    end
end
end

function L = optlag(x, y, tmax)
% Eq. 7.2: lag tau in [-tmax, tmax] maximizing corr(x(t), y(t+tau)).
% Lags are tried in the order 0, 1, -1, 2, -2, ... so ties go to the smallest |tau|.
best = -Inf;  L = 0;  n = numel(x);
for tau = [0, reshape([1:tmax; -(1:tmax)], 1, [])]
    if tau >= 0, a = x(1:n-tau);  b = y(1+tau:n);
    else,        a = x(1-tau:n);  b = y(1:n+tau); end
    c = corrcoef(a(:), b(:));  r = c(1,2);      % Pearson correlation (no toolbox needed)
    if r > best + 1e-12, best = r;  L = tau; end
end
end

function L = mean_abs_lag(agg, Id, members, tmax)
% Mean absolute optimal lag between an aggregate and its member counties.
L = mean(arrayfun(@(i) abs(optlag(agg, Id(:, i), tmax)), members));
end

function style_axes(ax, fs)
% Common axis style used in all figures.
set(ax, 'FontSize', fs, 'FontName', 'Times', 'LineWidth', 2, 'Box', 'on');
end

function M = make_row_stochastic(M)
% Reset each diagonal entry to 1 - (sum of the off-diagonal entries in its row),
% so that every row sums to one.
m = size(M, 1);
d = 1 - sum(M .* ~eye(m), 2);
if any(d < 0)
    warning('%d rows have off-diagonal sums above one; their diagonal is set to 0.', nnz(d < 0));
    d = max(d, 0);
end
M(logical(eye(m))) = d;
end

function save_fig(h, out_dir, name)
% Save a figure as .fig and .png in the result folders.
savefig(h, fullfile(out_dir, 'fig', [name '.fig']));
saveas(h,  fullfile(out_dir, 'png', [name '.png']));
end
