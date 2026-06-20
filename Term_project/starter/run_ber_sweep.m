% RUN_BER_SWEEP  BER vs Eb/N0 sweep harness, parameterised by Doppler.
%
% Produces one BER curve per Doppler value in DOPPLER_GRID. The default
% grid {0, 1000, 2000, 3000} Hz is the official term-project sweep:
%   - fd = 0    Hz : static multipath (Doppler diagnostic / sanity baseline)
%   - fd = 1000 Hz : moderate mobility
%   - fd = 2000 Hz : high mobility
%   - fd = 3000 Hz : the actual project operating point (FPV_C drone)
%
% Adjust N_TRIALS to trade speed vs precision while you iterate:
%   - Development (fast, noisy):   N_TRIALS = 20  – 50
%   - Final BER plot (smooth):     N_TRIALS = 200 – 500
%
% Each trial uses a different seed -> different impairment realisation.
% The grading harness is adaptive (>= 200 errors OR 5000 frames per point);
% you don't have to match that exactly, just make the final curve smooth.

clearvars; close all;

EbN0_dB       = 0:5:30;
DOPPLER_GRID  = [0, 1000, 2000, 3000];   % Hz
N_TRIALS      = 200;                     % final plot: smooth curves
seed0         = 1000;

% Tier-1 baseline at fd = 3000 Hz (from shipped stubs, measured)
tier1_ber = [4.0e-1, 2.2e-1, 1.3e-1, 1.0e-1, 8.5e-2, 7.9e-2, 7.53e-2];

BER = zeros(numel(DOPPLER_GRID), numel(EbN0_dB));

for d = 1:numel(DOPPLER_GRID)
    fd   = DOPPLER_GRID(d);
    opts = struct('fd_doppler', fd);
    fprintf('\n=== Doppler fd = %d Hz ===\n', fd);
    for i = 1:length(EbN0_dB)
        snr   = EbN0_dB(i);
        n_err = 0; n_bits = 0;
        for t = 1:N_TRIALS
            seed = seed0 + (d-1)*1e5 + (i-1)*N_TRIALS + t;
            [rx, params, ref] = ofdm_tx_and_channel(snr, seed, opts);
            bits = ofdm_rx(rx, params);
            L    = min(length(bits), length(ref));
            n_err  = n_err  + sum(bits(1:L) ~= ref(1:L));
            n_bits = n_bits + L;
        end
        BER(d, i) = n_err / max(n_bits, 1);
        fprintf('fd=%4d Hz  Eb/N0 = %4.1f dB    BER = %.4e   (%d errors / %d bits)\n', ...
                fd, snr, BER(d, i), n_err, n_bits);
    end
end

% ---- Plot ----------------------------------------------------------------
colors  = [0.22 0.45 0.70;   % blue   fd=0
           0.85 0.33 0.10;   % orange fd=1000
           0.47 0.67 0.19;   % green  fd=2000
           0.49 0.18 0.56];  % purple fd=3000
markers = {'o','s','^','d'};

figure('Position', [100 100 800 600]);
hold on;

for d = 1:numel(DOPPLER_GRID)
    semilogy(EbN0_dB, max(BER(d,:), 1e-6), ...
        ['-' markers{d}], ...
        'Color', colors(d,:), ...
        'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', colors(d,:), ...
        'DisplayName', sprintf('f_d = %d Hz', DOPPLER_GRID(d)));
end

% Tier-1 baseline fd=3000 Hz reference
semilogy(EbN0_dB, tier1_ber, '--k', ...
    'LineWidth', 1.5, 'MarkerSize', 7, ...
    'DisplayName', 'Tier-1 baseline  f_d = 3000 Hz');

set(gca, 'YScale', 'log');
grid on; grid minor;
xlabel('E_b/N_0 (dB)', 'FontSize', 12);
ylabel('BER',           'FontSize', 12);
title('Term Project — BER vs E_b/N_0 across Doppler grid', 'FontSize', 13);
ylim([1e-4 1]);
xlim([0 30]);
legend('Location', 'southwest', 'FontSize', 10);
set(gca, 'FontSize', 11);

% ---- Save as PNG with submission filename --------------------------------
LASTNAME = 'İlhan';   % <-- CHANGE THIS to your actual last name
out_png  = fullfile(fileparts(mfilename('fullpath')), '..', [LASTNAME '_ber.png']);
print(gcf, out_png, '-dpng', '-r150');
fprintf('\nSaved: %s\n', out_png);
