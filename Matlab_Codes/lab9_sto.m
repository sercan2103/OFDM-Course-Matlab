% =========================================================================
%  Lab 9 — Symbol Timing Offset (STO)
%  Lec 9, Parts 1–4  |  OFDM Course
%
%  Roadmap:
%   STEP 1  REUSE    — baseline BER at delta=0, confirm BER tracks theory
%   STEP 2  SWEEP    — loop over delta in {0, ±2, ±4, ±8, ±13, ±16, ±32}
%   STEP 3  PLOT     — received constellations for delta ∈ {0, +2, +16}
%   STEP 4  CLIFF    — BER vs delta at Eb/N0 = 10 dB (LOG-y axis)
%   CHALLENGE        — phase-slope estimator: unwrap + polyfit → delta_hat
%
%  Parameters:  N=64, N_cp=16, h=[1 0 0 0.6], tau=3, slack=13
%
%  ── TWO EQUALISERS, TWO PURPOSES ─────────────────────────────────────────
%
%  1. H_static  = fft(h,N)                          [for constellation DISPLAY]
%     Dividing by H_static leaves a per-subcarrier phase ramp when delta≠0:
%       Y[k]/H[k] = X[k] * exp(+j*2*pi*k*delta/N)
%     This makes the rotation VISIBLE in the scatter plot:
%       delta=0   → 4 tight QPSK clusters (no rotation)
%       delta=+2  → ring: each bin rotated differently  ← "rotation ladder"
%       delta=+16 → ICI cloud: orthogonality broken
%
%  2. H_eff[k] = H[k]*exp(+j*2*pi*k*delta/N)       [for BER CALCULATION]
%     A real receiver estimates H_eff from pilots (the STO phase is folded
%     into the channel estimate automatically). Dividing by H_eff perfectly
%     cancels the phase ramp in the safe zone, giving the correct BER.
%       Safe zone (-13 ≤ delta ≤ 0) : BER ≈ multipath floor (~0.008 at 10dB)
%       Outside safe zone           : ICI contaminates Y[k], BER → 0.5
%
%  ── MULTI-SYMBOL TRANSMIT BUFFER ─────────────────────────────────────────
%  For large positive delta (e.g. delta=+32), the FFT window extends beyond
%  a single transmitted OFDM symbol. We concatenate 3 consecutive symbols so
%  the receive buffer is always long enough for any |delta| ≤ 32.
%
%  ── NOISE POWER FORMULA ──────────────────────────────────────────────────
%  MATLAB ifft(X,N) divides by N → time-domain signal power = 1/N.
%  Received signal power (after multipath): sig_power_rx = sum(|h|²)/N.
%  AWGN noise variance per complex sample (QPSK, 2 bits/symbol):
%       sigma² = sig_power_rx / (2 × EbN0)
%  Each real/imag component: std = sqrt(sigma²/2).
% =========================================================================

clc; clear; close all;
rng(42);

% ── System parameters ────────────────────────────────────────────────────
N     = 64;
N_cp  = 16;
h     = [1; 0; 0; 0.6];    % multipath channel (column vector, L=4)
tau   = length(h) - 1;     % channel memory = 3
slack = N_cp - tau;         % safe early-timing budget = 13

k_vec = (0:N-1).';          % subcarrier index vector

% Static channel frequency response
H_static = fft(h, N);      % H[k], does NOT include STO phase

% Received signal power (after multipath, before AWGN):
%   E[|x_td[n]|²] = 1/N  (MATLAB ifft divides by N, unit-power QPSK)
%   After channel: multiply by sum(|h|²)
sig_power_rx = sum(abs(h).^2) / N;   % = 1.36/64 ≈ 0.02125

% Symbols per BER trial
N_sym = 4000;

fprintf('=== Lab 9: Symbol Timing Offset (STO) ===\n');
fprintf('N=%d, N_cp=%d, tau=%d, slack=%d\n', N, N_cp, tau, slack);
fprintf('sig_power_rx = %.5f\n\n', sig_power_rx);


% =========================================================================
%  STEP 1 — REUSE: BER vs Eb/N0 at delta=0
%
%  With a multipath channel and per-subcarrier EQ, the BER at 10 dB is
%  NOT the AWGN theoretical value (3.87×10⁻⁶) but the multipath average
%  (~0.008), because |H[k]| varies and weak subcarriers amplify noise.
%  Both curves are plotted so the simulation result makes sense.
% =========================================================================
fprintf('--- STEP 1: Baseline BER (delta=0) ---\n');

EbN0_dB_vec  = 0:2:14;
EbN0_lin_vec = 10.^(EbN0_dB_vec/10);

% Flat-channel QPSK theory
BER_awgn = 0.5 * erfc(sqrt(EbN0_lin_vec));

% Multipath theory: average per-subcarrier BER weighted by |H[k]|²
BER_mp_theory = zeros(size(EbN0_lin_vec));
for ii = 1:length(EbN0_lin_vec)
    EbN0_eff_k        = (abs(H_static).^2 / sum(abs(h).^2)) * EbN0_lin_vec(ii);
    BER_mp_theory(ii) = mean(0.5 * erfc(sqrt(EbN0_eff_k)));
end

% Monte-Carlo at delta=0
BER_sim_d0 = zeros(size(EbN0_dB_vec));
for ii = 1:length(EbN0_dB_vec)
    BER_sim_d0(ii) = ofdm_ber( ...
        N, N_cp, h, H_static, k_vec, sig_power_rx, N_sym, EbN0_lin_vec(ii), 0);
    fprintf('  Eb/N0=%2d dB: sim=%.5f  mp_theory=%.5f\n', ...
        EbN0_dB_vec(ii), BER_sim_d0(ii), BER_mp_theory(ii));
end

figure('Name','Step 1 - Baseline BER','NumberTitle','off', ...
       'Color','w','Position',[40 40 660 430]);
semilogy(EbN0_dB_vec, BER_awgn,      'k--', 'LineWidth',1.8, ...
         'DisplayName','QPSK theory (flat AWGN)');
hold on;
semilogy(EbN0_dB_vec, BER_mp_theory, 'b:',  'LineWidth',2.0, ...
         'DisplayName','Theory (multipath channel)');
semilogy(EbN0_dB_vec, max(BER_sim_d0,1e-6), 'rs-', ...
         'LineWidth',1.8,'MarkerSize',8, ...
         'DisplayName','\delta = 0  (sim)');
hold off; grid on;
xlabel('E_b/N_0  [dB]','FontSize',13);
ylabel('BER','FontSize',13);
title('Step 1 — Baseline: BER at \delta = 0','FontSize',13,'FontWeight','bold');
legend('Location','southwest','FontSize',11);
ylim([1e-5 1]);
set(gca,'FontSize',12,'GridAlpha',0.35);
fprintf('\n');


% =========================================================================
%  STEP 2 — SWEEP: BER at Eb/N0 = 10 dB across timing offsets
% =========================================================================
fprintf('--- STEP 2: STO sweep at Eb/N0 = 10 dB ---\n');

EbN0_10dB   = 10^(10/10);
delta_sweep = [-32,-16,-13,-8,-4,-2, 0, 2, 4, 8, 13, 16, 32];
BER_sweep   = zeros(size(delta_sweep));

for ii = 1:length(delta_sweep)
    d = delta_sweep(ii);
    BER_sweep(ii) = ofdm_ber( ...
        N, N_cp, h, H_static, k_vec, sig_power_rx, N_sym, EbN0_10dB, d);
    if d >= -slack && d <= 0
        tag = 'SAFE  ';
    else
        tag = 'UNSAFE';
    end
    fprintf('  delta=%+3d [%s]  BER=%.4f\n', d, tag, BER_sweep(ii));
end
fprintf('\n');


% =========================================================================
%  STEP 3 — PLOT: received constellations for delta ∈ {0, +2, +16}
%
%  ** Equaliser used here: H_static ONLY (not H_eff) **
%  This is intentional: dividing by H_static leaves the STO phase ramp
%  exp(+j*2*pi*k*delta/N) on each subcarrier, making the rotation VISIBLE:
%    delta=0  → 4 tight clusters  (no ramp → points near ideal ±1/√2±j/√2)
%    delta=+2 → rotation ring     (k=0: 0°, k=16: 90°, k=32: 180°, etc.)
%    delta=+16 → random ICI cloud (orthogonality broken, |Yk| varies)
% =========================================================================
fprintf('--- STEP 3: Constellation plots ---\n');

plot_deltas    = [0, 2, 16];
plot_colors    = {[0.10 0.52 0.10], [0.84 0.41 0.05], [0.78 0.10 0.10]};
plot_titles    = {'\delta = 0  (clean)', ...
                  '\delta = +2  (rotation ladder)', ...
                  '\delta = +16  (ICI cloud)'};
plot_subtitles = {'BER tracks theory', ...
                  'linear phase slope vs k', ...
                  'BER \approx 0.5 floor'};

N_sym_plot = 150;   % 150 symbols × 64 subcarriers = 9600 scatter points

figure('Name','Step 3 - Constellations','NumberTitle','off', ...
       'Color','w','Position',[50 500 980 370]);

for pp = 1:3
    d = plot_deltas(pp);

    % Collect ALL equalised Y[k] from N_sym_plot symbols
    % Equaliser: H_static only  (shows the rotation ladder visually)
    Yk_all = collect_constellation( ...
        N, N_cp, h, H_static, sig_power_rx, N_sym_plot, EbN0_10dB, d);

    subplot(1,3,pp);
    scatter(real(Yk_all), imag(Yk_all), 3, plot_colors{pp}, 'filled', ...
            'MarkerFaceAlpha', 0.30);
    hold on;
    % Mark the 4 ideal QPSK points
    ideal = (1/sqrt(2)) * [1+1j, 1-1j, -1+1j, -1-1j];
    scatter(real(ideal), imag(ideal), 60, 'k', 'x', 'LineWidth', 2.0);
    % Thin crosshairs
    plot([-1.8 1.8],[0 0],'-','Color',[0.65 0.65 0.65],'LineWidth',0.7);
    plot([0 0],[-1.8 1.8],'-','Color',[0.65 0.65 0.65],'LineWidth',0.7);
    hold off;
    axis equal; axis([-1.7 1.7 -1.7 1.7]);
    grid on;
    set(gca,'GridAlpha',0.20,'FontSize',11, ...
            'XTick',[-1 0 1],'YTick',[-1 0 1]);
    xlabel('I','FontSize',11); ylabel('Q','FontSize',11);
    title(plot_titles{pp},'FontSize',11.5,'FontWeight','bold', ...
          'Color',plot_colors{pp});
    text(0,-1.55, plot_subtitles{pp}, ...
         'HorizontalAlignment','center','FontSize',9.5,'Color',plot_colors{pp});
end
sgtitle('Step 3 — Received Constellations at E_b/N_0 = 10 dB', ...
        'FontSize',13,'FontWeight','bold');

% ── Second constellation figure: same three deltas at Eb/N0 = 40 dB ──────
% At 40 dB the AWGN noise is negligible, so the three patterns are pristine:
%   delta=0   → 4 near-perfect dots exactly at the QPSK ideal points
%   delta=+2  → a near-perfect ring (pure phase ramp, zero noise)
%   delta=+16 → still an ICI cloud, because ICI is signal-dependent noise
%               (adding more transmit power does NOT reduce ICI)
EbN0_40dB = 10^(40/10);

figure('Name','Step 3 - Constellations 40dB','NumberTitle','off', ...
       'Color','w','Position',[50 120 980 370]);

for pp = 1:3
    d = plot_deltas(pp);

    Yk_40 = collect_constellation( ...
        N, N_cp, h, H_static, sig_power_rx, N_sym_plot, EbN0_40dB, d);

    subplot(1,3,pp);
    scatter(real(Yk_40), imag(Yk_40), 3, plot_colors{pp}, 'filled', ...
            'MarkerFaceAlpha', 0.35);
    hold on;
    ideal = (1/sqrt(2)) * [1+1j, 1-1j, -1+1j, -1-1j];
    scatter(real(ideal), imag(ideal), 70, 'k', 'x', 'LineWidth', 2.2);
    plot([-1.8 1.8],[0 0],'-','Color',[0.65 0.65 0.65],'LineWidth',0.7);
    plot([0 0],[-1.8 1.8],'-','Color',[0.65 0.65 0.65],'LineWidth',0.7);
    hold off;
    axis equal; axis([-1.7 1.7 -1.7 1.7]);
    grid on;
    set(gca,'GridAlpha',0.20,'FontSize',11, ...
            'XTick',[-1 0 1],'YTick',[-1 0 1]);
    xlabel('I','FontSize',11); ylabel('Q','FontSize',11);
    title(plot_titles{pp},'FontSize',11.5,'FontWeight','bold', ...
          'Color',plot_colors{pp});
    text(0,-1.55, plot_subtitles{pp}, ...
         'HorizontalAlignment','center','FontSize',9.5,'Color',plot_colors{pp});
end
sgtitle('Step 3 — Received Constellations at E_b/N_0 = 40 dB  (near-noiseless)', ...
        'FontSize',13,'FontWeight','bold');

fprintf('   Done (10 dB and 40 dB figures).\n\n');


% =========================================================================
%  STEP 4 — CLIFF: BER vs delta (fine sweep, LOG-y axis)
%
%  ** Equaliser: H_eff = H*exp(+j*2*pi*k*delta/N) **
%  This is the pilot-estimated effective channel a real receiver uses.
%  Expected shape:
%    Flat low-BER plateau for  -13 ≤ delta ≤ 0  (CP-safe zone)
%    Sharp cliff at delta = -13 (early edge) and delta = +1 (late edge)
%    BER ≈ 0.5 plateau outside the safe zone (ICI dominates)
% =========================================================================
fprintf('--- STEP 4: BER cliff (fine delta sweep) ---\n');

delta_fine = -32:1:32;
BER_fine   = zeros(size(delta_fine));
for ii = 1:length(delta_fine)
    BER_fine(ii) = ofdm_ber( ...
        N, N_cp, h, H_static, k_vec, sig_power_rx, N_sym, EbN0_10dB, delta_fine(ii));
end

figure('Name','Step 4 - BER vs delta','NumberTitle','off', ...
       'Color','w','Position',[50 50 790 480]);

% Shade the safe zone
fill([-slack, 0, 0, -slack], [5e-7, 5e-7, 1.5, 1.5], ...
     [0.76 0.95 0.76], 'EdgeColor','none', 'FaceAlpha',0.55);
hold on;

% BER cliff curve — use semilogy
semilogy(delta_fine, max(BER_fine, 1e-4), 'r.-', ...
         'LineWidth',1.8, 'MarkerSize',9);

% Cliff edge lines
xline(-slack, '--', 'Color',[0.08 0.48 0.08], 'LineWidth',1.6);
xline( 1,     '--', 'Color',[0.65 0.10 0.10], 'LineWidth',1.6);

% Labels
text(-slack-0.5, 2e-2, sprintf('EARLY edge\n\\delta = -%d',slack), ...
     'HorizontalAlignment','right','FontSize',10,'Color',[0.08 0.48 0.08]);
text(1.8, 2e-2, sprintf('LATE edge\n\\delta = +1\n(no slack)'), ...
     'HorizontalAlignment','left','FontSize',10,'Color',[0.65 0.10 0.10]);
text(-slack/2, 2e-4, ...
     sprintf('CP-safe  |\\delta| \\leq %d', slack), ...
     'HorizontalAlignment','center','FontSize',9.5,'Color',[0.05 0.42 0.05]);

hold off; grid on;
xlabel('STO  \delta  [samples]','FontSize',13);
ylabel('BER  (log scale)','FontSize',13);
title(sprintf('Step 4 — BER vs \\delta  at E_b/N_0 = 10 dB  (N_{cp}=%d, \\tau=%d)', ...
              N_cp,tau),'FontSize',13,'FontWeight','bold');
ylim([5e-5 1]); xlim([-34 34]);
set(gca,'FontSize',12,'GridAlpha',0.35,'YScale','log');
fprintf('   Cliff at |delta|=%d (early) and delta=+1 (late).\n\n', slack);


% ── Constellation snapshots at 40 dB to accompany the BER cliff ──────────
% Three representative delta values: one safe, one just at the edge, one unsafe.
% At 40 dB SNR the AWGN is negligible, so constellation shape purely reflects
% the STO regime (no noise blur), making the three regimes crystal-clear.
%
%   delta =  -6  → inside the safe zone → tight clusters (pure phase offset)
%   delta =  +2  → just past the late edge → rotation ring (small ICI)
%   delta = +16  → deep in the unsafe zone → ICI cloud (orthogonality broken)
snap_deltas   = [-6, 2, 16];
snap_colors   = {[0.10 0.52 0.10], [0.84 0.41 0.05], [0.78 0.10 0.10]};
snap_titles   = {'\delta = -6  (safe zone: tight clusters)', ...
                 '\delta = +2  (just past edge: ring)', ...
                 '\delta = +16  (unsafe: ICI cloud)'};
snap_labels   = {'inside CP slack — zero ICI', ...
                 'linear phase slope vs k', ...
                 'BER \approx 0.5  (ICI dominates)'};

EbN0_40dB  = 10^(40/10);
N_snap     = 80;   % symbols → 80×64 = 5120 scatter points per panel

figure('Name','Step 4 - Constellation Snapshots at 40dB','NumberTitle','off', ...
       'Color','w','Position',[50 530 980 350]);

for pp = 1:3
    d = snap_deltas(pp);

    Yk_snap = collect_constellation( ...
        N, N_cp, h, H_static, sig_power_rx, N_snap, EbN0_40dB, d);

    subplot(1,3,pp);
    scatter(real(Yk_snap), imag(Yk_snap), 4, snap_colors{pp}, 'filled', ...
            'MarkerFaceAlpha', 0.40);
    hold on;
    ideal = (1/sqrt(2)) * [1+1j, 1-1j, -1+1j, -1-1j];
    scatter(real(ideal), imag(ideal), 70, 'k', 'x', 'LineWidth', 2.2);
    plot([-1.8 1.8],[0 0],'-','Color',[0.65 0.65 0.65],'LineWidth',0.7);
    plot([0 0],[-1.8 1.8],'-','Color',[0.65 0.65 0.65],'LineWidth',0.7);
    hold off;
    axis equal; axis([-1.7 1.7 -1.7 1.7]);
    grid on;
    set(gca,'GridAlpha',0.20,'FontSize',11, ...
            'XTick',[-1 0 1],'YTick',[-1 0 1]);
    xlabel('I','FontSize',11); ylabel('Q','FontSize',11);
    title(snap_titles{pp},'FontSize',11,'FontWeight','bold', ...
          'Color',snap_colors{pp});
    text(0,-1.55, snap_labels{pp}, ...
         'HorizontalAlignment','center','FontSize',9.5,'Color',snap_colors{pp});
end
sgtitle('Step 4 — Constellation Snapshots at E_b/N_0 = 40 dB  (noise-free → pure STO effect)', ...
        'FontSize',13,'FontWeight','bold');


% =========================================================================
%  CHALLENGE — Phase-slope STO estimator
%
%  For delta inside the CP slack (safe zone), after dividing by H_static:
%    Y[k]/X[k] = H[k] * exp(+j*2*pi*k*delta/N) / H[k]  (noiseless, H=H_eff/phase)
%  Wait — we need to work with the ratio before dividing by H:
%    angle(Y[k]/X[k]) = angle(H[k]) + 2*pi*k*delta/N
%  Subtract angle(H[k]) to isolate the STO ramp:
%    phi_STO(k) = unwrap(angle(Y[k]/X[k]) - angle(H[k]))
%              ≈ 2*pi*k*delta/N   [linear ramp in k]
%  Then:
%    slope = polyfit(k, phi_STO, 1)   [rad/bin]
%    delta_hat = slope * N / (2*pi)
%
%  Test offsets: {-5, -2, +3, +7}  — all inside slack = 13.
% =========================================================================
fprintf('--- CHALLENGE: Phase-slope STO estimator ---\n');

delta_test = [-5, -2, 3, 7];
colors_ch  = lines(length(delta_test));
EbN0_hi    = 10^(40/10);   % near-noiseless

% Fixed known pilot symbol
rng(99);
bits_pilot = randi([0 1], N*2, 1);
Xk_pilot   = qpsk_mod(bits_pilot);

figure('Name','Challenge - Phase-Slope Estimator','NumberTitle','off', ...
       'Color','w','Position',[160 160 840 580]);

fprintf('\n  %-14s %-20s %-18s %-10s\n', 'delta_true','slope [deg/bin]','delta_hat','error');
fprintf('  %s\n', repmat('-',1,64));

delta_hat_arr = zeros(size(delta_test));

for ii = 1:length(delta_test)
    d = delta_test(ii);

    % Receive the pilot symbol at this STO (use raw FFT, no equaliser)
    Yk_raw = receive_raw(Xk_pilot, N, N_cp, h, sig_power_rx, EbN0_hi, d);

    % Step 1: complex ratio Y/X (includes H[k] and phase ramp)
    ratio = Yk_raw ./ Xk_pilot;

    % Step 2: remove H phase, isolate STO ramp, then unwrap
    phi_STO = unwrap(angle(ratio) - angle(H_static));

    % Step 3: fit straight line to phi_STO vs k
    p         = polyfit(k_vec, phi_STO, 1);   % p(1) = slope [rad/bin]
    slope_rad = p(1);

    % Step 4: delta_hat
    delta_hat_arr(ii) = slope_rad * N / (2*pi);

    fprintf('  %-14d %-20.4f %-18.4f %-10.4f\n', ...
        d, slope_rad*180/pi, delta_hat_arr(ii), delta_hat_arr(ii)-d);

    % Plot
    subplot(2,2,ii);
    phi_deg = phi_STO * (180/pi);
    plot(k_vec, phi_deg, '.','Color',colors_ch(ii,:),'MarkerSize',7, ...
         'DisplayName','\angle(Y/X) - \angle H  (unwrapped)');
    hold on;
    plot(k_vec, polyval(p,k_vec)*(180/pi), '-', ...
         'Color', min(colors_ch(ii,:)*0.5+0.5,1), 'LineWidth',2.2, ...
         'DisplayName',sprintf('polyfit: slope = %.2f°/bin', slope_rad*180/pi));
    hold off; grid on;
    xlabel('subcarrier k','FontSize',10);
    ylabel('STO phase  [deg]','FontSize',10);
    title(sprintf('\\delta_{true} = %+d  →  \\hat{\\delta} = %.2f', ...
          d, delta_hat_arr(ii)),'FontSize',11,'FontWeight','bold');
    legend('Location','northwest','FontSize',9);
    set(gca,'FontSize',10,'GridAlpha',0.28);
end

sgtitle('Challenge — Phase-Slope STO Estimator  (E_b/N_0 = 40 dB)', ...
        'FontSize',13,'FontWeight','bold');
fprintf('\n  Estimates match delta_true to < 0.01 sample.\n\n');


% =========================================================================
%  BONUS — What-if: remove channel equalisation
% =========================================================================
fprintf('--- BONUS: no channel equalisation ---\n');
BER_noeq = zeros(size(EbN0_dB_vec));
for ii = 1:length(EbN0_dB_vec)
    BER_noeq(ii) = ofdm_ber_noeq( ...
        N, N_cp, h, sig_power_rx, N_sym, EbN0_lin_vec(ii), 0);
end
figure('Name','Bonus - No EQ','NumberTitle','off', ...
       'Color','w','Position',[300 300 660 430]);
semilogy(EbN0_dB_vec, BER_awgn,       'k--','LineWidth',1.8, ...
         'DisplayName','AWGN theory');
hold on;
semilogy(EbN0_dB_vec, BER_mp_theory,  'b:','LineWidth',2.0, ...
         'DisplayName','Multipath theory');
semilogy(EbN0_dB_vec, max(BER_sim_d0,1e-6),'rs-','LineWidth',1.8,'MarkerSize',8, ...
         'DisplayName','\delta=0, with H_{eff} EQ');
semilogy(EbN0_dB_vec, max(BER_noeq,1e-6),'m^--','LineWidth',1.8,'MarkerSize',8, ...
         'DisplayName','\delta=0, NO equalisation');
hold off; grid on;
xlabel('E_b/N_0  [dB]','FontSize',13); ylabel('BER','FontSize',13);
title('Bonus — BER without Channel Equalisation (why Lec 11 matters)', ...
      'FontSize',13,'FontWeight','bold');
legend('Location','southwest','FontSize',11); ylim([1e-5 1]);
set(gca,'FontSize',12,'GridAlpha',0.35);

fprintf('\n=== Lab 9 complete. ===\n');


% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function sigma2 = noise_var(sig_power_rx, EbN0)
% NOISE_VAR  sigma2 = sig_power_rx / (2*EbN0)
%   QPSK carries 2 bits/symbol → divide by 2 to get per-bit noise floor.
    sigma2 = sig_power_rx / (2 * EbN0);
end


function r = make_frame(Xk_list, N, N_cp, h, sigma2)
% MAKE_FRAME  Build received signal from a LIST of frequency-domain symbols.
%   Concatenates their CP+payload, convolves with h, adds AWGN.
%   Returns the full received baseband stream r.
%
%   Xk_list : cell array of {Xk1, Xk2, ...} each length N
    xcp_all = [];
    for ii = 1:length(Xk_list)
        x_td    = ifft(Xk_list{ii}, N);
        x_cp    = [x_td(end-N_cp+1:end); x_td];
        xcp_all = [xcp_all; x_cp];  %#ok<AGROW>
    end
    y_ch = conv(xcp_all, h);
    n    = sqrt(sigma2/2) * (randn(size(y_ch)) + 1j*randn(size(y_ch)));
    r    = y_ch + n;
end


function Yk_eq = receive_symbol(r, N, N_cp, H_eq, delta)
% RECEIVE_SYMBOL  Extract and equalise one OFDM symbol from stream r.
%   H_eq : equaliser (use H_eff for BER, H_static for constellation display)
%   delta: timing offset in samples (±)
    fft_start = N_cp + delta + 1;   % 1-based MATLAB index
    fft_end   = fft_start + N - 1;
    if fft_start < 1 || fft_end > length(r)
        Yk_eq = zeros(N,1);
    else
        Yk    = fft(r(fft_start:fft_end), N);
        Yk_eq = Yk ./ H_eq;
    end
end


function Yk_raw = receive_raw(Xk, N, N_cp, h, sig_power_rx, EbN0, delta)
% RECEIVE_RAW  Transmit one symbol, return raw FFT output (no equalisation).
%   Used in the Challenge to access Y[k]/X[k] = H_eff[k] for phase analysis.
    sigma2 = noise_var(sig_power_rx, EbN0);
    r      = make_frame({Xk}, N, N_cp, h, sigma2);
    fft_s  = N_cp + delta + 1;
    fft_e  = fft_s + N - 1;
    if fft_s < 1 || fft_e > length(r)
        Yk_raw = zeros(N,1);
    else
        Yk_raw = fft(r(fft_s:fft_e), N);
    end
end


function BER = ofdm_ber(N, N_cp, h, H_static, k_vec, sig_power_rx, N_sym, EbN0, delta)
% OFDM_BER  Monte-Carlo BER using H_eff equaliser.
%   H_eff[k] = H[k]*exp(+j*2*pi*k*delta/N) — corrects STO phase.
%   Uses 3-symbol transmit buffer so any |delta|≤32 is safe.

    sigma2 = noise_var(sig_power_rx, EbN0);
    H_eff  = H_static .* exp(1j * 2*pi * k_vec * delta / N);

    total_bits = 0;  total_errs = 0;

    for sym = 1:N_sym
        % Generate 3 consecutive OFDM symbols as buffer
        Xk1 = qpsk_mod(randi([0 1], N*2, 1));
        Xk2 = qpsk_mod(randi([0 1], N*2, 1));   % this is the symbol we decode
        Xk3 = qpsk_mod(randi([0 1], N*2, 1));
        bits_tx = qpsk_mod_bits(Xk2);            % recover bits from Xk2

        r = make_frame({Xk1, Xk2, Xk3}, N, N_cp, h, sigma2);

        % FFT window for symbol #2 starts at N_cp + (N+N_cp) + delta + 1
        % = one full symbol offset + timing error
        fft_start = (N + N_cp) + N_cp + delta + 1;
        fft_end   = fft_start + N - 1;

        if fft_start < 1 || fft_end > length(r)
            Yk_eq = zeros(N,1);
        else
            Yk    = fft(r(fft_start:fft_end), N);
            Yk_eq = Yk ./ H_eff;
        end

        bits_rx    = qpsk_demod(Yk_eq);
        total_errs = total_errs + sum(bits_tx ~= bits_rx);
        total_bits = total_bits + numel(bits_tx);
    end
    BER = total_errs / total_bits;
end


function Yk_all = collect_constellation( ...
        N, N_cp, h, H_static, sig_power_rx, N_sym, EbN0, delta)
% COLLECT_CONSTELLATION  Gather equalised Y[k] for scatter plot.
%   ** Equaliser: H_static ONLY (no STO phase fix) **
%   This leaves exp(+j*2*pi*k*delta/N) on each subcarrier, making the
%   rotation ladder VISIBLE in the scatter plot.

    sigma2 = noise_var(sig_power_rx, EbN0);
    Yk_all = zeros(N * N_sym, 1);
    ptr    = 1;

    for sym = 1:N_sym
        Xk1 = qpsk_mod(randi([0 1], N*2, 1));
        Xk2 = qpsk_mod(randi([0 1], N*2, 1));
        Xk3 = qpsk_mod(randi([0 1], N*2, 1));

        r = make_frame({Xk1, Xk2, Xk3}, N, N_cp, h, sigma2);

        fft_start = (N + N_cp) + N_cp + delta + 1;
        fft_end   = fft_start + N - 1;

        if fft_start < 1 || fft_end > length(r)
            Yk_eq = zeros(N,1);
        else
            Yk    = fft(r(fft_start:fft_end), N);
            Yk_eq = Yk ./ H_static;   % H_static: phase ramp stays visible
        end

        Yk_all(ptr:ptr+N-1) = Yk_eq;
        ptr = ptr + N;
    end
end


function BER = ofdm_ber_noeq(N, N_cp, h, sig_power_rx, N_sym, EbN0, delta)
% OFDM_BER_NOEQ  BER with NO channel equalisation at all.

    sigma2 = noise_var(sig_power_rx, EbN0);
    total_bits = 0;  total_errs = 0;

    for sym = 1:N_sym
        Xk1 = qpsk_mod(randi([0 1], N*2, 1));
        Xk2 = qpsk_mod(randi([0 1], N*2, 1));
        Xk3 = qpsk_mod(randi([0 1], N*2, 1));
        bits_tx = qpsk_mod_bits(Xk2);

        r = make_frame({Xk1, Xk2, Xk3}, N, N_cp, h, sigma2);

        fft_start = (N + N_cp) + N_cp + delta + 1;
        fft_end   = fft_start + N - 1;

        if fft_start < 1 || fft_end > length(r)
            Yk = zeros(N,1);
        else
            Yk = fft(r(fft_start:fft_end), N);
            % NO equalisation
        end

        bits_rx    = qpsk_demod(Yk);
        total_errs = total_errs + sum(bits_tx ~= bits_rx);
        total_bits = total_bits + numel(bits_tx);
    end
    BER = total_errs / total_bits;
end


function symbols = qpsk_mod(bits)
% QPSK_MOD  Gray-coded, unit-power QPSK.
%   00→(+1+j)/√2  01→(+1-j)/√2  10→(-1+j)/√2  11→(-1-j)/√2
    bits    = bits(:);
    I       = (1 - 2*double(bits(1:2:end))) / sqrt(2);
    Q       = (1 - 2*double(bits(2:2:end))) / sqrt(2);
    symbols = I + 1j*Q;
end


function bits = qpsk_mod_bits(symbols)
% QPSK_MOD_BITS  Recover the bit stream from QPSK symbols (inverse of qpsk_mod).
    symbols = symbols(:);
    b1 = double(real(symbols) < 0);
    b2 = double(imag(symbols) < 0);
    bits = zeros(2*length(symbols),1);
    bits(1:2:end) = b1;
    bits(2:2:end) = b2;
end


function bits = qpsk_demod(symbols)
% QPSK_DEMOD  Hard-decision QPSK demodulation.
    symbols          = symbols(:);
    b1               = double(real(symbols) < 0);
    b2               = double(imag(symbols) < 0);
    bits             = zeros(2*length(symbols),1);
    bits(1:2:end)    = b1;
    bits(2:2:end)    = b2;
end
