% =========================================================================
%  ofdm_M4_equalizer.m
%  M4 Method: Simple ICI-Aware Equalizer for OFDM under CFO + Doppler
%
%  System: N=64, N_cp=16, QPSK, h=[1 0 0 0.6], epsilon=0.05, f_d~240 Hz
%
%  Produces
%    Figure 1 — BER vs Eb/N0: no correction | M4 | theoretical AWGN
%    Figure 2 — Constellation before / after M4 at Eb/N0 = 15 dB
%
%  =========================================================================
%  KEY DESIGN DECISIONS (why previous versions failed)
%  =========================================================================
%
%  1. NOISE CALIBRATION
%     MATLAB's ifft uses the 1/N convention, so the time-domain signal
%     power = (1/N) * frequency-domain power.  Adding AWGN of variance
%     sigma2 per sample, then taking the N-point FFT, gives noise power
%     N*sigma2 per subcarrier (FFT sums N samples coherently for signal
%     but incoherently for noise -> noise power scales by N).
%     Signal power per subcarrier = 1 (normalised QPSK, freq domain).
%
%     Correct Es/N0 relationship:
%       Es/N0 = 1 / (N * sigma2)  =>  sigma2 = 1 / (N * bps * Eb/N0)
%
%     Previous error: sigma2 = 1/(2*bps*Eb/N0)  [missed the 1/N factor]
%     Effect: noise was 32x too strong, burying all signal.
%
%  2. CFO + DOPPLER MODEL
%     Physical sample rate: fs = N * Delta_f = 64 * 15e3 = 960 kHz
%     CFO:     per-sample phase = 2*pi*epsilon / N
%     Doppler: per-sample phase = 2*pi*f_d / fs = 2*pi*f_d / (N*Delta_f)
%     Combined eps_eff = epsilon + f_d/Delta_f = 0.05 + 0.016 = 0.066
%     Combined per-sample phase dphi = 2*pi*eps_eff / N
%
%  3. CPE REMOVAL (Step 1 of M4)
%     After the N-point FFT, the exact expression is:
%       Y[k] = exp(j*dphi*(n0+N_cp)) * sum_m X[m]*H[m]*sinc(m+eps_eff-k)
%                                             *exp(j*pi*(m+eps_eff-k)*(N-1)/N)
%     The CPE is exp(j*dphi*(n0+N_cp)) — the phase at the FIRST useful
%     sample (NOT the centre of the window).
%     Previous error: used (n0 + N_cp + (N-1)/2) — added an extra
%     pi*eps_eff*(N-1)/N rotation that then conflicted with Step 2.
%
%  4. CHANNEL + SINC EQUALISATION (Step 2 of M4)
%     After CPE removal, the self-term on subcarrier k is:
%       X[k]*H[k]*sinc(eps_eff)*exp(j*pi*eps_eff*(N-1)/N)
%     Divide by H[k]*sinc(eps_eff)*exp(j*pi*eps_eff*(N-1)/N) = H[k]*c_full
%     to recover X[k] (plus residual ICI).
%     c_full = sinc(eps_eff) * exp(j*pi*eps_eff*(N-1)/N)
%
%  5. ICI CANCELLATION (Step 3 of M4)
%     Dominant ICI on k comes from m=k±1 with coefficient:
%       alpha = sinc(eps_eff-1) / sinc(eps_eff)   (ratio normalised by c_self)
%     Hard-decision estimates Xhat from Step 2 give the interferers.
% =========================================================================

clear; clc; close all;

%% ── System Parameters ────────────────────────────────────────────────────
N        = 64;           % Subcarriers
N_cp     = 16;           % Cyclic prefix
N_sym    = 1000;         % OFDM symbols per Eb/N0 point
bps      = 2;            % Bits per symbol (QPSK)

% Physical parameters
Delta_f  = 15e3;                        % Subcarrier spacing [Hz]
fs       = N * Delta_f;                 % Sample rate = 960 kHz
f_c      = 2.4e9;                       % Carrier [Hz]
v_ms     = 30;                          % Speed [m/s]
c_light  = 3e8;
f_d      = (v_ms / c_light) * f_c;     % Doppler shift ~ 240 Hz

% Normalised offsets
epsilon  = 0.05;                        % CFO (normalised to Delta_f)
eps_dop  = f_d / Delta_f;              % Doppler normalised = 0.016
eps_eff  = epsilon + eps_dop;           % Combined = 0.066

% Per-sample phase increment
dphi     = 2 * pi * eps_eff / N;        % rad/sample

% Channel
h_cir    = [1, 0, 0, 0.6];

% Eb/N0 sweep
EbN0_dB       = [0, 5, 10, 15, 20, 25];
EbN0_dB_const = 15;                     % For constellation plot

fprintf('=== M4 ICI-Aware OFDM Equalizer Simulation ===\n');
fprintf('N=%d, N_cp=%d, fs=%.0f kHz, Delta_f=%.0f Hz\n', N, N_cp, fs/1e3, Delta_f);
fprintf('epsilon=%.4f, f_d=%.1f Hz, eps_dop=%.4f, eps_eff=%.4f\n', ...
        epsilon, f_d, eps_dop, eps_eff);
fprintf('dphi = 2*pi*eps_eff/N = %.6f rad/sample\n', dphi);
fprintf('Channel h = ['); fprintf('%.2f ', h_cir); fprintf(']\n');

%% ── QPSK Mapping ─────────────────────────────────────────────────────────
qpsk_map = (1/sqrt(2)) * [1+1j, -1+1j, 1-1j, -1-1j];

%% ── Frequency-Domain Channel ─────────────────────────────────────────────
H_freq = fft(h_cir, N);

%% ── M4 Equaliser Coefficients ────────────────────────────────────────────
% c_full: the exact complex attenuation of the desired subcarrier after
%         CPE removal (self-term of the ICI sum)
c_full = sinc(eps_eff) * exp(1j * pi * eps_eff * (N-1) / N);

% ICI cancellation coefficient (nearest-neighbour leakage ratio)
alpha  = sinc(eps_eff - 1) / sinc(eps_eff);

fprintf('\nM4 coefficients:\n');
fprintf('  c_full = sinc(%.4f)*exp(j*pi*%.4f*(N-1)/N) = %.6f angle %.2f deg\n', ...
        eps_eff, eps_eff, abs(c_full), angle(c_full)*180/pi);
fprintf('  alpha  = sinc(%.4f-1)/sinc(%.4f) = %.6f\n\n', eps_eff, eps_eff, alpha);

%% ── Allocate Results ─────────────────────────────────────────────────────
BER_no_corr  = zeros(1, length(EbN0_dB));
BER_M4       = zeros(1, length(EbN0_dB));
BER_theory   = zeros(1, length(EbN0_dB));
rx_before_eq = [];    % constellation storage (after CPE only)
rx_after_eq  = [];    % constellation storage (after full M4)

%% ═══════════════════════════════════════════════════════════════════════
%%  MAIN SIMULATION LOOP
%% ═══════════════════════════════════════════════════════════════════════
for snr_idx = 1:length(EbN0_dB)

    EbN0_lin = 10^(EbN0_dB(snr_idx) / 10);

    % ── Noise Variance ──────────────────────────────────────────────────
    % Time-domain signal power = (1/N) per sample (MATLAB ifft convention).
    % After N-point FFT, signal power = 1 per subcarrier.
    % Noise added in time domain (power sigma2/sample) → after FFT: N*sigma2/subcarrier.
    % For correct Eb/N0: 1 / (N * sigma2) = bps * Eb/N0
    %   => sigma2 = 1 / (N * bps * Eb/N0)
    sigma2 = 1 / (N * bps * EbN0_lin);

    % Theoretical QPSK BER (AWGN only)
    BER_theory(snr_idx) = qfunc(sqrt(2 * EbN0_lin));

    n_err_no = 0;
    n_err_M4 = 0;
    n_bits   = 0;

    for sym_idx = 0:(N_sym - 1)

        %% ── Transmitter ─────────────────────────────────────────────────
        bits_tx = randi([0, 3], 1, N);
        X       = qpsk_map(bits_tx + 1);
        x_td    = ifft(X, N);                        % time-domain (power 1/N per sample)
        x_cp    = [x_td(end-N_cp+1:end), x_td];     % add cyclic prefix

        %% ── Multipath Channel ────────────────────────────────────────────
        x_mp = conv(x_cp, h_cir);
        x_mp = x_mp(1 : N + N_cp);

        %% ── CFO + Doppler (per-sample exponential) ───────────────────────
        % Absolute sample indices for this symbol's block
        n0    = sym_idx * (N + N_cp);
        n_vec = n0 : n0 + (N + N_cp) - 1;
        r_imp = x_mp .* exp(1j * dphi .* n_vec);

        %% ── AWGN ─────────────────────────────────────────────────────────
        % Complex noise: real and imag each ~ N(0, sigma2)
        noise = sqrt(sigma2) * (randn(1, N+N_cp) + 1j*randn(1, N+N_cp));
        r_rx  = r_imp + noise;

        %% ── CP Removal + FFT ─────────────────────────────────────────────
        r_useful = r_rx(N_cp + 1 : end);    % strip CP
        Y        = fft(r_useful, N);         % N-point FFT

        %% ── Baseline Decoder (no correction) ────────────────────────────
        [~, idx_no]  = min(abs(repmat(Y.', 1, 4) - repmat(qpsk_map, N, 1)), [], 2);
        bits_rx_no   = idx_no' - 1;
        n_err_no     = n_err_no + sum(bits_tx ~= bits_rx_no);

        %% ─────────────────────────────────────────────────────────────────
        %%  M4 EQUALIZER
        %% ─────────────────────────────────────────────────────────────────

        %% Step 1 — CPE Removal
        % After FFT, Y[k] contains exp(j*dphi*(n0+N_cp)) as a common factor.
        % Remove it: Y1 = Y * exp(-j*dphi*(n0+N_cp))
        % NOTE: use (n0+N_cp) = first useful sample index, NOT the centre.
        phi_CPE = dphi * (n0 + N_cp);
        Y1      = Y * exp(-1j * phi_CPE);

        %% Step 2 — Channel + Sinc Equalisation
        % Self-term of Y1[k] = X[k]*H[k]*c_full
        % Divide by H[k]*c_full to recover X[k] (+ICI residual +noise)
        Y2 = Y1 ./ (H_freq * c_full);

        %% Step 3 — Nearest-Neighbour ICI Cancellation
        % Hard-decision estimates of X from Step 2
        [~, idx_tmp] = min(abs(repmat(Y2.', 1, 4) - repmat(qpsk_map, N, 1)), [], 2);
        Xhat         = qpsk_map(idx_tmp');

        % ICI from subcarrier k-1 and k+1 (circular)
        Xhat_m1 = circshift(Xhat,  1);    % Xhat[k-1]
        Xhat_p1 = circshift(Xhat, -1);    % Xhat[k+1]

        % Subtract estimated ICI contributions
        Y3 = Y2 - alpha * (Xhat_m1 + Xhat_p1);

        %% M4 Final Decision
        [~, idx_M4] = min(abs(repmat(Y3.', 1, 4) - repmat(qpsk_map, N, 1)), [], 2);
        bits_rx_M4  = idx_M4' - 1;
        n_err_M4    = n_err_M4 + sum(bits_tx ~= bits_rx_M4);

        n_bits = n_bits + N;

        %% Store constellation at target Eb/N0
        if EbN0_dB(snr_idx) == EbN0_dB_const
            rx_before_eq = [rx_before_eq, Y1];   %#ok
            rx_after_eq  = [rx_after_eq,  Y3];   %#ok
        end

    end % symbol loop

    BER_no_corr(snr_idx) = n_err_no / n_bits;
    BER_M4(snr_idx)      = n_err_M4 / n_bits;

    fprintf('Eb/N0=%2d dB | no-corr=%7.4f | M4=%7.4f | theory=%.6f\n', ...
        EbN0_dB(snr_idx), BER_no_corr(snr_idx), BER_M4(snr_idx), BER_theory(snr_idx));

end % Eb/N0 loop

%% ═══════════════════════════════════════════════════════════════════════
%%  FIGURE 1 — BER vs Eb/N0
%% ═══════════════════════════════════════════════════════════════════════
figure(1); clf;
semilogy(EbN0_dB, BER_no_corr, 'r-o',  'LineWidth', 2, 'MarkerSize', 8, ...
         'DisplayName', 'No correction (baseline)');
hold on;
semilogy(EbN0_dB, BER_M4,      'b-s',  'LineWidth', 2, 'MarkerSize', 8, ...
         'DisplayName', 'M4 ICI-aware equalizer');
semilogy(EbN0_dB, BER_theory,  'k--^', 'LineWidth', 2, 'MarkerSize', 8, ...
         'DisplayName', 'Theoretical AWGN (QPSK)');
hold off;
grid on;
xlabel('Eb/N0 (dB)', 'FontSize', 13);
ylabel('BER',        'FontSize', 13);
title({'OFDM BER vs Eb/N0 — CFO + Doppler Impairments', ...
       sprintf('N=%d, N_{cp}=%d, \\epsilon=%.2f, f_d=%.0f Hz', ...
               N, N_cp, epsilon, f_d)}, 'FontSize', 13);
legend('Location', 'southwest', 'FontSize', 11);
ylim([1e-4, 1]);
set(gca, 'FontSize', 11);
annotation('textbox', [0.14, 0.14, 0.42, 0.14], ...
    'String', {sprintf('h=[1 0 0 0.6], \\epsilon=%.2f, f_d=%.0f Hz', epsilon, f_d), ...
               sprintf('\\epsilon_{eff} = \\epsilon + f_d/\\Deltaf = %.4f', eps_eff), ...
               sprintf('M4: CPE + H[k]\\cdotc_{full} + ICI cancel (\\alpha=%.4f)', alpha)}, ...
    'FitBoxToText', 'on', 'BackgroundColor', [0.97 0.97 0.97], ...
    'EdgeColor', [0.5 0.5 0.5], 'FontSize', 9);

%% ═══════════════════════════════════════════════════════════════════════
%%  FIGURE 2 — Constellation at Eb/N0 = 15 dB
%% ═══════════════════════════════════════════════════════════════════════
figure(2); clf;
ideal_pts = qpsk_map;

subplot(1, 2, 1);
scatter(real(rx_before_eq), imag(rx_before_eq), 3, [0.35 0.35 0.85], 'filled', ...
        'MarkerFaceAlpha', 0.15);
hold on;
scatter(real(ideal_pts), imag(ideal_pts), 200, 'r', 'x', 'LineWidth', 2.5);
hold off;
grid on; axis equal;
xlabel('In-Phase'); ylabel('Quadrature');
title({'Before M4 Equalizer', ...
       sprintf('(after CPE only, Eb/N_0=%d dB)', EbN0_dB_const)}, 'FontSize', 11);
legend('Received', 'Ideal QPSK', 'Location', 'northeast', 'FontSize', 9);
xlim([-2.5, 2.5]); ylim([-2.5, 2.5]);
set(gca, 'FontSize', 10);

subplot(1, 2, 2);
scatter(real(rx_after_eq), imag(rx_after_eq), 3, [0.05 0.50 0.10], 'filled', ...
        'MarkerFaceAlpha', 0.15);
hold on;
scatter(real(ideal_pts), imag(ideal_pts), 200, 'r', 'x', 'LineWidth', 2.5);
hold off;
grid on; axis equal;
xlabel('In-Phase'); ylabel('Quadrature');
title({'After M4 Equalizer', ...
       sprintf('(CPE + H[k] + ICI cancel, Eb/N_0=%d dB)', EbN0_dB_const)}, 'FontSize', 11);
legend('Received', 'Ideal QPSK', 'Location', 'northeast', 'FontSize', 9);
xlim([-2.5, 2.5]); ylim([-2.5, 2.5]);
set(gca, 'FontSize', 10);

sgtitle({'QPSK Constellation: Before vs After M4 Equalizer', ...
         sprintf('\\epsilon=%.2f, f_d=%.0f Hz, h=[1 0 0 0.6]', epsilon, f_d)}, ...
        'FontSize', 13);

fprintf('\nDone. Figures 1 and 2 ready.\n');
