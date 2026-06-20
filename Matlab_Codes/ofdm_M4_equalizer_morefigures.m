% =========================================================================
%  ofdm_M4_equalizer.m
%  M4 Method: Simple ICI-Aware Equalizer for OFDM under CFO + Doppler
%
%  System: N=64, N_cp=16, QPSK, h=[1 0 0 0.6], epsilon=0.05, f_d~240 Hz
%
%  FIGURES PRODUCED
%  ----------------
%  Figure 1  — BER vs Eb/N0 (three curves)
%  Figure 2  — Constellation before / after M4  (Eb/N0 = 15 dB)
%  Figure 3  — M4 pipeline: 4-step constellation progression
%  Figure 4  — ICI leakage spectrum (sinc model)
%  Figure 5  — ICI coefficient matrix (16×16 heat-map)
%  Figure 6  — Channel frequency response |H[k]| and angle
%  Figure 7  — CPE phase drift vs symbol index
%  Figure 8  — Residual interference power at each M4 step (bar chart)
%  Figure 9  — EVM² vs Eb/N0 at each M4 step (waterfall)
%  Figure 10 — ICI coefficient alpha vs eps_eff (sensitivity)
%
% =========================================================================

clear; clc; close all;

% ── System Parameters ──────────────────────────────────────────────────
N        = 64;
N_cp     = 16;
N_sym    = 1000;
bps      = 2;

Delta_f  = 15e3;
fs       = N * Delta_f;
f_c      = 2.4e9;
v_ms     = 30;
c_light  = 3e8;
f_d      = (v_ms / c_light) * f_c;

epsilon  = 0.05;
eps_dop  = f_d / Delta_f;
eps_eff  = epsilon + eps_dop;
dphi     = 2 * pi * eps_eff / N;

h_cir    = [1, 0, 0, 0.6];

EbN0_dB       = [0, 5, 10, 15, 20, 25];
EbN0_dB_const = 15;

fprintf('=== M4 ICI-Aware OFDM Equalizer Simulation ===\n');
fprintf('N=%d, N_cp=%d, fs=%.0f kHz\n', N, N_cp, fs/1e3);
fprintf('epsilon=%.4f, f_d=%.1f Hz, eps_eff=%.4f\n', epsilon, f_d, eps_eff);

% ── Derived quantities ─────────────────────────────────────────────────
qpsk_map = (1/sqrt(2)) * [1+1j, -1+1j, 1-1j, -1-1j];
H_freq   = fft(h_cir, N);
c_full   = sinc(eps_eff) * exp(1j * pi * eps_eff * (N-1) / N);
alpha    = sinc(eps_eff - 1) / sinc(eps_eff);

fprintf('c_full = %.4f angle %.2f deg\n', abs(c_full), angle(c_full)*180/pi);
fprintf('alpha  = %.6f\n\n', alpha);

% ── Allocate BER arrays ────────────────────────────────────────────────
BER_no_corr = zeros(1, length(EbN0_dB));
BER_M4      = zeros(1, length(EbN0_dB));
BER_theory  = zeros(1, length(EbN0_dB));

% ── Extra storage for diagnostic plots ────────────────────────────────
rx_raw_all   = [];   % Y  (no correction)       — for Fig 3
rx_Y1_all    = [];   % Y1 (after CPE)           — for Fig 3
rx_Y2_all    = [];   % Y2 (after channel EQ)    — for Fig 3
rx_Y3_all    = [];   % Y3 (after ICI cancel)    — for Fig 3

% Residual power accumulators (no-noise run, done after main loop)
% EVM vs SNR storage
EVM_Y1_dB   = zeros(1, length(EbN0_dB));
EVM_Y2_dB   = zeros(1, length(EbN0_dB));
EVM_Y3_dB   = zeros(1, length(EbN0_dB));

% ═══════════════════════════════════════════════════════════════════════
%  MAIN SIMULATION LOOP
% ═══════════════════════════════════════════════════════════════════════
for snr_idx = 1:length(EbN0_dB)

    EbN0_lin = 10^(EbN0_dB(snr_idx) / 10);
    sigma2   = 1 / (N * bps * EbN0_lin);
    BER_theory(snr_idx) = qfunc(sqrt(2 * EbN0_lin));

    n_err_no = 0;  n_err_M4 = 0;  n_bits = 0;
    evm1 = 0;  evm2 = 0;  evm3 = 0;   % MSE accumulators for EVM

    for sym_idx = 0:(N_sym - 1)

        % TX
        bits_tx = randi([0, 3], 1, N);
        X       = qpsk_map(bits_tx + 1);
        x_td    = ifft(X, N);
        x_cp    = [x_td(end-N_cp+1:end), x_td];

        % Multipath
        x_mp = conv(x_cp, h_cir);
        x_mp = x_mp(1 : N + N_cp);

        % CFO + Doppler
        n0    = sym_idx * (N + N_cp);
        n_vec = n0 : n0 + (N + N_cp) - 1;
        r_imp = x_mp .* exp(1j * dphi .* n_vec);

        % AWGN
        noise = sqrt(sigma2) * (randn(1,N+N_cp) + 1j*randn(1,N+N_cp));
        r_rx  = r_imp + noise;

        % CP strip + FFT
        Y = fft(r_rx(N_cp+1:end), N);

        % Baseline
        [~, idx_no] = min(abs(repmat(Y.',1,4) - repmat(qpsk_map,N,1)), [], 2);
        n_err_no    = n_err_no + sum(bits_tx ~= idx_no'-1);

        % ── M4 Equalizer ──────────────────────────────────────────────

        % Step 1: CPE removal
        phi_CPE = dphi * (n0 + N_cp);
        Y1      = Y * exp(-1j * phi_CPE);

        % Step 2: Channel + sinc equalisation
        Y2 = Y1 ./ (H_freq * c_full);

        % Step 3: Nearest-neighbour ICI cancellation
        [~, idx_tmp] = min(abs(repmat(Y2.',1,4) - repmat(qpsk_map,N,1)), [], 2);
        Xhat    = qpsk_map(idx_tmp');
        Y3      = Y2 - alpha * (circshift(Xhat,1) + circshift(Xhat,-1));

        % Final decision
        [~, idx_M4] = min(abs(repmat(Y3.',1,4) - repmat(qpsk_map,N,1)), [], 2);
        n_err_M4    = n_err_M4 + sum(bits_tx ~= idx_M4'-1);
        n_bits      = n_bits + N;

        % EVM accumulation (MSE relative to true X)
        evm1 = evm1 + sum(abs(Y1 - X .* H_freq * c_full).^2);
        evm2 = evm2 + sum(abs(Y2 - X).^2);
        evm3 = evm3 + sum(abs(Y3 - X).^2);

        % Store constellations at target Eb/N0 (first 200 symbols)
        if EbN0_dB(snr_idx) == EbN0_dB_const && sym_idx < 200
            rx_raw_all = [rx_raw_all, Y ];  %#ok
            rx_Y1_all  = [rx_Y1_all,  Y1]; %#ok
            rx_Y2_all  = [rx_Y2_all,  Y2]; %#ok
            rx_Y3_all  = [rx_Y3_all,  Y3]; %#ok
        end

    end % symbol loop

    BER_no_corr(snr_idx) = n_err_no / n_bits;
    BER_M4(snr_idx)      = n_err_M4 / n_bits;
    EVM_Y1_dB(snr_idx)   = 10*log10(evm1 / n_bits);
    EVM_Y2_dB(snr_idx)   = 10*log10(evm2 / n_bits);
    EVM_Y3_dB(snr_idx)   = 10*log10(evm3 / n_bits);

    fprintf('Eb/N0=%2d dB | no-corr=%.4f | M4=%.4f | theory=%.6f\n', ...
        EbN0_dB(snr_idx), BER_no_corr(snr_idx), BER_M4(snr_idx), BER_theory(snr_idx));

end % Eb/N0 loop

% ── No-noise run to measure pure ICI residual at each step ────────────
step_pow_raw = 0;  step_pow_Y2 = 0;  step_pow_Y3 = 0;
N_nn = 500;
for sym_idx = 0:(N_nn-1)
    bits_tx = randi([0,3],1,N);
    X       = qpsk_map(bits_tx+1);
    x_td    = ifft(X,N);
    x_cp    = [x_td(end-N_cp+1:end), x_td];
    x_mp    = conv(x_cp, h_cir); x_mp = x_mp(1:N+N_cp);
    n0      = sym_idx*(N+N_cp);
    r       = x_mp .* exp(1j*dphi*(n0:n0+N+N_cp-1));
    Y       = fft(r(N_cp+1:end), N);
    Y1      = Y * exp(-1j * dphi*(n0+N_cp));
    Y2      = Y1 ./ (H_freq * c_full);
    [~,itmp]= min(abs(repmat(Y2.',1,4)-repmat(qpsk_map,N,1)),[],2);
    Xh      = qpsk_map(itmp');
    Y3      = Y2 - alpha*(circshift(Xh,1)+circshift(Xh,-1));
    step_pow_raw = step_pow_raw + mean(abs(Y/N - X).^2);  % normalised
    step_pow_Y2  = step_pow_Y2  + mean(abs(Y2 - X).^2);
    step_pow_Y3  = step_pow_Y3  + mean(abs(Y3 - X).^2);
end
step_pow_raw = step_pow_raw / N_nn;
step_pow_Y2  = step_pow_Y2  / N_nn;
step_pow_Y3  = step_pow_Y3  / N_nn;

fprintf('\nNo-noise ICI residual power:\n');
fprintf('  Raw (no EQ):    %.4f  (%.1f dB)\n', step_pow_raw, 10*log10(step_pow_raw));
fprintf('  After Step 2:   %.4f  (%.1f dB)\n', step_pow_Y2,  10*log10(step_pow_Y2));
fprintf('  After Step 3:   %.4f  (%.1f dB)\n', step_pow_Y3,  10*log10(step_pow_Y3));

fprintf('\nAll figures being generated...\n');

% ═══════════════════════════════════════════════════════════════════════
%  FIGURE 1 — BER vs Eb/N0
% ═══════════════════════════════════════════════════════════════════════
figure(1); clf;
semilogy(EbN0_dB, BER_no_corr, 'r-o',  'LineWidth',2,'MarkerSize',8, ...
         'DisplayName','No correction (baseline)');
hold on;
semilogy(EbN0_dB, BER_M4,      'b-s',  'LineWidth',2,'MarkerSize',8, ...
         'DisplayName','M4 ICI-aware equalizer');
semilogy(EbN0_dB, BER_theory,  'k--^', 'LineWidth',2,'MarkerSize',8, ...
         'DisplayName','Theoretical AWGN (QPSK)');
hold off;
grid on;
xlabel('Eb/N0 (dB)','FontSize',13); ylabel('BER','FontSize',13);
title({'OFDM BER vs Eb/N0 — CFO + Doppler Impairments', ...
       sprintf('N=%d,  N_{cp}=%d,  \\epsilon=%.2f,  f_d=%.0f Hz', ...
               N,N_cp,epsilon,f_d)},'FontSize',13);
legend('Location','southwest','FontSize',11);
ylim([1e-4,1]);
annotation('textbox',[0.14,0.13,0.44,0.14], ...
    'String',{sprintf('h=[1 0 0 0.6],  \\epsilon=%.2f,  f_d=%.0f Hz',epsilon,f_d), ...
              sprintf('\\epsilon_{eff} = \\epsilon + f_d/\\Deltaf = %.4f',eps_eff), ...
              sprintf('M4: CPE + H[k]\\cdotc_{full} + ICI cancel  (\\alpha=%.4f)',alpha)}, ...
    'FitBoxToText','on','BackgroundColor',[0.97 0.97 0.97], ...
    'EdgeColor',[0.5 0.5 0.5],'FontSize',9);

% ═══════════════════════════════════════════════════════════════════════
%  FIGURE 2 — Constellation before / after M4  (original required plot)
% ═══════════════════════════════════════════════════════════════════════
figure(2); clf;
subplot(1,2,1);
scatter(real(rx_Y1_all),imag(rx_Y1_all),3,[0.35 0.35 0.85],'filled', ...
        'MarkerFaceAlpha',0.15);
hold on;
scatter(real(qpsk_map),imag(qpsk_map),200,'r','x','LineWidth',2.5);
hold off; grid on; axis equal;
xlabel('In-Phase'); ylabel('Quadrature');
title({'Before M4 Equalizer', ...
       sprintf('(after CPE only,  Eb/N_0=%d dB)',EbN0_dB_const)},'FontSize',11);
legend('Received','Ideal QPSK','Location','northeast','FontSize',9);
xlim([-2.5 2.5]); ylim([-2.5 2.5]); set(gca,'FontSize',10);

subplot(1,2,2);
scatter(real(rx_Y3_all),imag(rx_Y3_all),3,[0.05 0.50 0.10],'filled', ...
        'MarkerFaceAlpha',0.15);
hold on;
scatter(real(qpsk_map),imag(qpsk_map),200,'r','x','LineWidth',2.5);
hold off; grid on; axis equal;
xlabel('In-Phase'); ylabel('Quadrature');
title({'After M4 Equalizer', ...
       sprintf('(CPE + H[k] + ICI cancel,  Eb/N_0=%d dB)',EbN0_dB_const)},'FontSize',11);
legend('Received','Ideal QPSK','Location','northeast','FontSize',9);
xlim([-2.5 2.5]); ylim([-2.5 2.5]); set(gca,'FontSize',10);
sgtitle({'QPSK Constellation: Before vs After M4 Equalizer', ...
         sprintf('\\epsilon=%.2f,  f_d=%.0f Hz,  h=[1 0 0 0.6]',epsilon,f_d)}, ...
        'FontSize',13);

% ═══════════════════════════════════════════════════════════════════════
%  FIGURE 3 — M4 Pipeline: 4-stage constellation progression
% ═══════════════════════════════════════════════════════════════════════
figure(3); clf;
stage_data  = {rx_raw_all, rx_Y1_all, rx_Y2_all, rx_Y3_all};
stage_title = {'Step 0: Raw FFT output', ...
               'Step 1: After CPE removal', ...
               'Step 2: After channel EQ', ...
               'Step 3: After ICI cancel (M4)'};
stage_note  = {'Uniform cloud — ICI + phase drift', ...
               'Phase stabilised, channel still in', ...
               'Clusters centred, ICI residual spread', ...
               'Tight clusters \rightarrow reliable decisions'};
stage_col   = {[0.45 0.45 0.45]; [0.20 0.35 0.80]; [0.85 0.45 0.10]; [0.10 0.60 0.20]};
lim_vals    = {[-4 4]; [-4 4]; [-2 2]; [-2 2]};

for s = 1:4
    subplot(1,4,s);
    rx = stage_data{s};
    scatter(real(rx), imag(rx), 2, stage_col{s}, 'filled', 'MarkerFaceAlpha', 0.12);
    hold on;
    scatter(real(qpsk_map), imag(qpsk_map), 120, 'r', 'x', 'LineWidth', 2.2);
    hold off; grid on; axis equal;
    xlim(lim_vals{s}); ylim(lim_vals{s});
    title(stage_title{s}, 'FontSize', 9.5, 'FontWeight', 'bold');
    xlabel('In-Phase', 'FontSize', 8);
    if s == 1; ylabel('Quadrature', 'FontSize', 8); end
    text(0.5, -0.10, stage_note{s}, 'Units','normalized', ...
         'HorizontalAlignment','center','FontSize',7.5, ...
         'Color', stage_col{s}, 'FontAngle','italic');
    % decision region boundaries (dashed diagonals)
    line(lim_vals{s}, [0 0], 'Color',[0.7 0.7 0.7],'LineStyle',':','LineWidth',0.8);
    line([0 0], lim_vals{s}, 'Color',[0.7 0.7 0.7],'LineStyle',':','LineWidth',0.8);
end
sgtitle(sprintf('M4 Equalizer Pipeline — Constellation at Each Step  (Eb/N_0 = %d dB)', ...
                EbN0_dB_const), 'FontSize', 12, 'FontWeight', 'bold');

% ═══════════════════════════════════════════════════════════════════════
%  FIGURE 4 — ICI Leakage Spectrum (sinc model)
% ═══════════════════════════════════════════════════════════════════════
figure(4); clf;
d_arr  = linspace(-4, 4, 800);
sinc_v = abs(sinc(eps_eff - d_arr));   % |sinc(eps_eff - Delta_k)|

plot(d_arr, sinc_v, 'b-', 'LineWidth', 2); hold on;
% Shade self-term region
fill([eps_eff-0.3, eps_eff+0.3, eps_eff+0.3, eps_eff-0.3], [0,0,1.05,1.05], ...
     [0.2 0.6 0.9], 'FaceAlpha',0.20,'EdgeColor','none');
% Shade nearest neighbours
fill([eps_eff-1.3, eps_eff-0.7, eps_eff-0.7, eps_eff-1.3],[0,0,1.05,1.05], ...
     [1 0.5 0.1],'FaceAlpha',0.25,'EdgeColor','none');
fill([eps_eff+0.7, eps_eff+1.3, eps_eff+1.3, eps_eff+0.7],[0,0,1.05,1.05], ...
     [1 0.5 0.1],'FaceAlpha',0.25,'EdgeColor','none');

xline(eps_eff,   '--', 'Color',[0.2 0.6 0.2], 'LineWidth',1.5, ...
      'Label', sprintf(' Self (k=m, \\Delta k=\\epsilon_{eff}=%.3f)',eps_eff), ...
      'LabelHorizontalAlignment','right','FontSize',8);
xline(eps_eff-1, '--', 'Color',[0.85 0.4 0.0], 'LineWidth',1.5, ...
      'Label',' k-1  ICI', 'LabelHorizontalAlignment','right','FontSize',8);
xline(eps_eff+1, '--', 'Color',[0.85 0.4 0.0], 'LineWidth',1.5, ...
      'Label',' k+1  ICI', 'LabelHorizontalAlignment','left','FontSize',8);

% Annotate values
text(eps_eff+0.08, sinc(0)+0.02, ...
     sprintf('|sinc(%.3f)| = %.3f', eps_eff, sinc(eps_eff)), ...
     'FontSize',8,'Color',[0.1 0.5 0.1]);
text(eps_eff-1+0.08, abs(sinc(eps_eff-1))+0.04, ...
     sprintf('\\alpha = %.4f', abs(alpha)), ...
     'FontSize',8,'Color',[0.7 0.3 0.0]);

hold off; grid on;
xlabel('\Delta k  (subcarrier offset from target)', 'FontSize', 12);
ylabel('|sinc(\epsilon_{eff} - \Delta k)|',          'FontSize', 12);
title({'ICI Leakage Spectrum — sinc Model', ...
       sprintf('\\epsilon_{eff}=%.4f,  nearest-neighbour \\alpha=%.4f', eps_eff, abs(alpha))}, ...
      'FontSize', 12);
xlim([-4 4]); ylim([0 1.1]);
legend('sinc leakage envelope','Self region','k±1 ICI (M4 cancels)', ...
       'Location','north','FontSize',9);

% ═══════════════════════════════════════════════════════════════════════
%  FIGURE 5 — ICI Coefficient Matrix (16×16 heat-map)
% ═══════════════════════════════════════════════════════════════════════
figure(5); clf;
Nshow = 16;
ICI_mat = zeros(Nshow, Nshow);
for k = 1:Nshow
    for m = 1:Nshow
        delta_km = (m-1) + eps_eff - (k-1);
        % wrap into (-N/2, N/2]
        delta_km = mod(delta_km + N/2, N) - N/2;
        ICI_mat(k,m) = abs(sinc(delta_km) * exp(1j*pi*delta_km*(N-1)/N) / c_full);
    end
end
imagesc(0:Nshow-1, 0:Nshow-1, ICI_mat);
colormap(flipud(gray));  colorbar;
hold on;
% Overlay diagonal (self) in cyan and ±1 in orange
for k = 0:Nshow-1
    plot(k, k, 'cs', 'MarkerSize', 9, 'MarkerFaceColor','c', 'LineWidth',0.5);
    if k > 0
        plot(k-1, k, 'o', 'Color',[1 0.45 0.0], 'MarkerSize',6, ...
             'MarkerFaceColor',[1 0.75 0.3], 'LineWidth',0.5);
    end
    if k < Nshow-1
        plot(k+1, k, 'o', 'Color',[1 0.45 0.0], 'MarkerSize',6, ...
             'MarkerFaceColor',[1 0.75 0.3], 'LineWidth',0.5);
    end
end
hold off;
xlabel('Source subcarrier  m', 'FontSize', 12);
ylabel('Target subcarrier  k', 'FontSize', 12);
title({'ICI Coefficient Matrix  |c(k,m)|  (first 16 subcarriers)', ...
       'Diagonal = self-term (≈1),  Off-diagonal = leakage M4 cancels'}, 'FontSize',11);
legend({'Self  k=m', 'k±1  (M4 cancels)'}, 'Location','eastoutside','FontSize',9);
set(gca,'XTick',0:2:Nshow-1,'YTick',0:2:Nshow-1);
cb = colorbar; cb.Label.String = '|coefficient|'; cb.Label.FontSize = 9;

% ═══════════════════════════════════════════════════════════════════════
%  FIGURE 6 — Channel Frequency Response
% ═══════════════════════════════════════════════════════════════════════
figure(6); clf;
k_arr = 0:N-1;

subplot(2,1,1);
bar(k_arr, abs(H_freq), 0.75, 'FaceColor',[0.22 0.45 0.75], 'EdgeColor','none');
hold on;
yline(1,'--','Color',[0.5 0.5 0.5],'LineWidth',1.2,'Label','Flat reference','FontSize',8);
% Shade fades
idx_fade = find(abs(H_freq) < 0.8);
for ii = idx_fade
    patch([ii-1.5 ii-0.5 ii-0.5 ii-1.5],[0 0 2 2],[1 0.3 0.3], ...
          'FaceAlpha',0.15,'EdgeColor','none');
end
hold off; grid on;
ylabel('|H[k]|','FontSize',11);
title(sprintf('Channel Frequency Response  h = [1, 0, 0, 0.6]  (N=%d)',N),'FontSize',12);
xlim([-1 N]); ylim([0 1.85]);
text(1,1.70,sprintf('max=%.2f,  min=%.2f,  causes %.1f dB variation', ...
     max(abs(H_freq)), min(abs(H_freq)), ...
     20*log10(max(abs(H_freq))/min(abs(H_freq)))),'FontSize',8,'Color',[0.4 0.4 0.4]);

subplot(2,1,2);
plot(k_arr, unwrap(angle(H_freq))*180/pi, '-', ...
     'Color',[0.55 0.15 0.75],'LineWidth',1.8);
grid on;
xlabel('Subcarrier index  k','FontSize',11);
ylabel('\angleH[k]  (degrees)','FontSize',11);
title('Channel Phase Response (unwrapped)','FontSize',11);
xlim([-1 N]);

% ═══════════════════════════════════════════════════════════════════════
%  FIGURE 7 — CPE Phase Drift vs Symbol Index
% ═══════════════════════════════════════════════════════════════════════
figure(7); clf;
sym_show    = 0:99;
phi_drift   = mod(dphi * (sym_show*(N+N_cp) + N_cp) * 180/pi, 360);
phi_cumul   = dphi * (sym_show*(N+N_cp) + N_cp) * 180/pi;

subplot(2,1,1);
plot(sym_show, phi_drift, 'Color',[0.85 0.35 0.10], 'LineWidth',1.8);
hold on;
% QPSK boundary lines at 45,135,225,315
for bnd = [45 135 225 315]
    yline(bnd,'--','Color',[0.8 0.2 0.2],'LineWidth',1,'Alpha',0.6);
end
scatter(sym_show, phi_drift, 18, 'filled', 'MarkerFaceColor',[0.85 0.35 0.10]);
hold off; grid on;
xlabel('OFDM Symbol Index','FontSize',11);
ylabel('CPE Phase mod 360°  (deg)','FontSize',11);
title({'Common Phase Error (CPE) Across Symbols — Before Step 1 Correction', ...
       sprintf('\\epsilon_{eff}=%.4f,  phase increment = %.4f deg/symbol', ...
               eps_eff, dphi*(N+N_cp)*180/pi)},'FontSize',11);
ylim([0 370]);
text(2,330,{'Dashed lines = QPSK decision boundaries (±45°,±135°)', ...
            'Phase crosses boundaries \Rightarrow errors without correction'}, ...
    'FontSize',8,'Color',[0.7 0.1 0.1]);

subplot(2,1,2);
plot(sym_show, phi_cumul, 'Color',[0.20 0.40 0.80], 'LineWidth',1.8);
hold on;
scatter(sym_show, phi_cumul, 18, 'filled', 'MarkerFaceColor',[0.20 0.40 0.80]);
hold off; grid on;
xlabel('OFDM Symbol Index','FontSize',11);
ylabel('Cumulative CPE (degrees)','FontSize',11);
title('Cumulative Phase Growth (unbounded without correction)','FontSize',11);

% ═══════════════════════════════════════════════════════════════════════
%  FIGURE 8 — Residual Interference Power at Each M4 Step (bar chart)
% ═══════════════════════════════════════════════════════════════════════
figure(8); clf;
step_labels = {'No EQ (raw)', 'After Step 1\n(CPE only)', ...
               'After Step 2\n(+ Channel EQ)', 'After Step 3\n(+ ICI cancel)'};
step_powers_dB = [10*log10(step_pow_raw + 1), ...   % +1 avoids log(0) for raw
                  10*log10(step_pow_raw + 1), ...   % CPE alone doesn't change ICI
                  10*log10(step_pow_Y2), ...
                  10*log10(step_pow_Y3)];
% Recompute raw and CPE-only residuals cleanly
raw_pow  = step_pow_raw;   % measured above
cpe_pow  = step_pow_raw;   % CPE is a scalar, doesn't change ICI magnitude
vals_dB  = [10*log10(raw_pow), 10*log10(cpe_pow), ...
            10*log10(step_pow_Y2), 10*log10(step_pow_Y3)];
bar_colors = [0.45 0.45 0.45; 0.80 0.35 0.35; 0.90 0.55 0.10; 0.15 0.60 0.25];

for b = 1:4
    bar_h = bar(b, vals_dB(b), 0.55);
    bar_h.FaceColor = bar_colors(b,:);
    hold on;
end
hold off; grid on;
set(gca, 'XTick', 1:4, 'XTickLabel', ...
    {'Raw FFT','Step 1: CPE','Step 2: +Channel','Step 3: +ICI'}, 'FontSize',10);
ylabel('Residual Interference Power (dB)','FontSize',12);
title({'Residual ICI Power at Each M4 Stage  (no noise, ICI only)', ...
       'Lower is better — Step 2 provides the largest gain'},'FontSize',12);
for b = 1:4
    text(b, vals_dB(b)+0.4, sprintf('%.1f dB', vals_dB(b)), ...
         'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
end
yline(0,'--','Color',[0.5 0.5 0.5],'LineWidth',1,'Label','0 dB ref','FontSize',8);

% ═══════════════════════════════════════════════════════════════════════
%  FIGURE 9 — EVM² vs Eb/N0 at each M4 step  (waterfall)
% ═══════════════════════════════════════════════════════════════════════
figure(9); clf;
plot(EbN0_dB, EVM_Y1_dB, 'r-o', 'LineWidth',2,'MarkerSize',8, ...
     'DisplayName','After Step 1: CPE removal');
hold on;
plot(EbN0_dB, EVM_Y2_dB, '-s', 'Color',[0.90 0.50 0.05], 'LineWidth',2,'MarkerSize',8, ...
     'DisplayName','After Step 2: + Channel EQ');
plot(EbN0_dB, EVM_Y3_dB, 'g-^', 'LineWidth',2,'MarkerSize',8, ...
     'DisplayName','After Step 3: + ICI cancel  (M4 output)');

% Shade improvement bands
fill([EbN0_dB, fliplr(EbN0_dB)], ...
     [EVM_Y1_dB, fliplr(EVM_Y2_dB)], ...
     [0.95 0.75 0.55], 'FaceAlpha',0.25,'EdgeColor','none', ...
     'DisplayName','Gain from Step 2');
fill([EbN0_dB, fliplr(EbN0_dB)], ...
     [EVM_Y2_dB, fliplr(EVM_Y3_dB)], ...
     [0.65 0.90 0.65], 'FaceAlpha',0.25,'EdgeColor','none', ...
     'DisplayName','Gain from Step 3');
hold off; grid on;

xlabel('Eb/N0 (dB)',   'FontSize',12);
ylabel('Residual Error Power  E[|Y - X|^2]  (dB)','FontSize',12);
title({'Error Vector Magnitude² at Each M4 Step vs Eb/N0', ...
       'Shows which step dominates the correction at each SNR'},'FontSize',12);
legend('Location','northeast','FontSize',9);

% ═══════════════════════════════════════════════════════════════════════
%  FIGURE 10 — ICI coefficient alpha vs eps_eff  (sensitivity)
% ═══════════════════════════════════════════════════════════════════════
figure(10); clf;
eps_sweep = linspace(0.001, 0.49, 500);
alpha_sweep = abs(sinc(eps_sweep - 1) ./ sinc(eps_sweep));
sinc_att_dB = 20*log10(sinc(eps_sweep));    % sinc self-attenuation

yyaxis left;
plot(eps_sweep, alpha_sweep, 'Color',[0.55 0.15 0.75], 'LineWidth',2.2);
ylabel('\alpha  =  |sinc(\epsilon-1) / sinc(\epsilon)|', ...
       'Color',[0.55 0.15 0.75],'FontSize',11);
ylim([0 0.65]);

hold on;
% Mark operating point
xline(eps_eff,'--','Color',[0.85 0.40 0.05],'LineWidth',1.8, ...
      'Label',sprintf(' \\epsilon_{eff}=%.3f\\newline\\alpha=%.3f', eps_eff, abs(alpha)), ...
      'LabelHorizontalAlignment','right','FontSize',9);
scatter(eps_eff, abs(alpha), 80, [0.85 0.40 0.05], 'filled', 'ZData',1);
hold off;

yyaxis right;
plot(eps_sweep, sinc_att_dB, 'Color',[0.15 0.55 0.30], 'LineWidth',1.8, 'LineStyle','--');
ylabel('sinc self-attenuation  (dB)', 'Color',[0.15 0.55 0.30],'FontSize',11);
ylim([-5 0.2]);

grid on;
xlabel('\epsilon_{eff}  (normalised frequency offset)','FontSize',12);
title({'M4 Sensitivity: ICI Coefficient \alpha and sinc Attenuation vs \epsilon_{eff}', ...
       'Both corrections embedded in c_{full} and \alpha — computed once, fixed for the run'}, ...
      'FontSize',11);
legend({'\alpha (left axis,  ICI cancel strength)', ...
        'sinc attenuation (right axis, Step 2 compensates)'}, ...
       'Location','northwest','FontSize',9);

fprintf('\nAll 10 figures generated.\n');
