clear; close all; clc;
rng(6);

%% Parameters
N_BITS   = 5000;
SPS      = 8;
ALPHA    = 0.35;
N_TAPS   = 65;
EbN0_dB  = 0:2:12;
BITS_PER_SYM = 1; % BPSK

rrc = get_rrc_filter(SPS, ALPHA, N_TAPS);
gd = (numel(rrc) - 1) / 2;

%% Common Tx chain
bits = randi([0 1], 1, N_BITS);
symbols = 2*bits - 1; % BPSK mapping: 0->-1, 1->+1

tx_up = zeros(1, N_BITS*SPS);
tx_up(1:SPS:end) = symbols;
tx = conv(tx_up, rrc, 'full');

%% TODO 1 - AWGN Tx/Rx chain (with matched filter)
ber_awgn = zeros(size(EbN0_dB));
ber_theory = zeros(size(EbN0_dB));

%% TODO 2 - 2-tap multipath channel
h1 = [1, 0,0, 0.9]; % echo at 1T delay (8 samples when SPS=8)
h2 = [1, 0,0, 0.6]; % echo at 1T delay (8 samples when SPS=8)
tx_mp = conv(tx, h1, 'full');
tx_mp2 = conv(tx, h2, 'full');
ber_h1 = zeros(size(EbN0_dB));
ber_h2 = zeros(size(EbN0_dB));

% Store one operating point for eye/MF plots
EYE_EBN0_DB = 8;
idx_eye = find(EbN0_dB == EYE_EBN0_DB, 1);
eyeSig_awgn_mf = [];
eyeSig_h1_mf = [];
eyeSig_h2_mf = [];
mfSamp_awgn = [];
mfSamp_h1 = [];
mfSamp_h2 = [];

for ii = 1:numel(EbN0_dB)
	ebn0 = EbN0_dB(ii);

	% AWGN-only channel
	rx_awgn = add_awgn_from_ebn0(tx, ebn0, SPS, BITS_PER_SYM);
	y_awgn_mf = conv(rx_awgn, rrc, 'full');

	start_awgn = 2*gd + 1;
	z_awgn = y_awgn_mf(start_awgn:SPS:start_awgn + (N_BITS-1)*SPS);
	b_hat_awgn = z_awgn > 0;
	ber_awgn(ii) = mean(b_hat_awgn ~= bits);

	gamma_b = 10^(ebn0/10);
	ber_theory(ii) = 0.5 * erfc(sqrt(gamma_b));

	% 2-tap multipath + AWGN
	rx_mp = add_awgn_from_ebn0(tx_mp, ebn0, SPS, BITS_PER_SYM);
	y_mp_mf = conv(rx_mp, rrc, 'full');

	rx_mp2 = add_awgn_from_ebn0(tx_mp2, ebn0, SPS, BITS_PER_SYM);
	y_mp2_mf = conv(rx_mp2, rrc, 'full');

	% Main path is at 0 delay; echo creates ISI but does not shift decision origin.
	start_mp = 2*gd + 1;
	z_mp = y_mp_mf(start_mp:SPS:start_mp + (N_BITS-1)*SPS);
	b_hat_mp = z_mp > 0;
	ber_h1(ii) = mean(b_hat_mp ~= bits);

	z_mp2 = y_mp2_mf(start_mp:SPS:start_mp + (N_BITS-1)*SPS);
	b_hat_mp2 = z_mp2 > 0;
	ber_h2(ii) = mean(b_hat_mp2 ~= bits);

	if ~isempty(idx_eye) && ii == idx_eye
		eyeSig_awgn_mf = y_awgn_mf;
		eyeSig_h1_mf = y_mp_mf;
		eyeSig_h2_mf = y_mp2_mf;
		mfSamp_awgn = z_awgn;
		mfSamp_h1 = z_mp;
		mfSamp_h2 = z_mp2;
	end

	fprintf('Eb/N0=%2d dB | BER(AWGN)=%.3e | Theory=%.3e | BER(h1)=%.3e | BER(h2)=%.3e\n', ...
		ebn0, ber_awgn(ii), ber_theory(ii), ber_h1(ii), ber_h2(ii));
end

%% TODO 3 - Comparison plots
figure('Name', 'Lab 6 BER Comparison', 'Color', 'w');
semilogy(EbN0_dB, ber_theory, 'k--', 'LineWidth', 1.8); hold on;
semilogy(EbN0_dB, ber_awgn, 'o-', 'LineWidth', 2.0, 'MarkerSize', 7, 'Color', [0.10 0.45 0.75]);
semilogy(EbN0_dB, ber_h1, 's-', 'LineWidth', 2.0, 'MarkerSize', 7, 'Color', [0.85 0.33 0.10]);
semilogy(EbN0_dB, ber_h2, 'd-', 'LineWidth', 2.0, 'MarkerSize', 7, 'Color', [0.20 0.65 0.20]);
grid on;
xlabel('Eb/N0 (dB)');
ylabel('BER');
title('BER: AWGN vs 2-tap Multipath (h1 and h2)');
legend('BPSK Theory', 'Simulation (AWGN + MF)', 'Simulation (h1 + MF)', 'Simulation (h2 + MF)', 'Location', 'southwest');
ylim([1e-5 1]);

if ~isempty(eyeSig_awgn_mf) && ~isempty(eyeSig_h1_mf) && ~isempty(eyeSig_h2_mf)
	figure('Name', 'Lab 6 Eye Diagram Comparison', 'Color', 'w');
	subplot(1,3,1);
	plot_eye_custom(eyeSig_awgn_mf, SPS, 220, 2*gd + 1, [0.10 0.45 0.75]);
	grid on; xlim([0 2]);
	xlabel('Time (T)'); ylabel('Amplitude');
	title(sprintf('Eye @ MF Output (AWGN), Eb/N0 = %d dB', EYE_EBN0_DB));

	subplot(1,3,2);
	plot_eye_custom(eyeSig_h1_mf, SPS, 220, 2*gd + 1, [0.85 0.33 0.10]);
	grid on; xlim([0 2]);
	xlabel('Time (T)'); ylabel('Amplitude');
	title(sprintf('Eye @ MF Output (h1), Eb/N0 = %d dB', EYE_EBN0_DB));

	subplot(1,3,3);
	plot_eye_custom(eyeSig_h2_mf, SPS, 220, 2*gd + 1, [0.20 0.65 0.20]);
	grid on; xlim([0 2]);
	xlabel('Time (T)'); ylabel('Amplitude');
	title(sprintf('Eye @ MF Output (h2), Eb/N0 = %d dB', EYE_EBN0_DB));
end

if ~isempty(mfSamp_awgn) && ~isempty(mfSamp_h1) && ~isempty(mfSamp_h2)
	nShow = min(120, N_BITS);
	k = 1:nShow;
	figure('Name', 'Lab 6 Matched-Filter Sample Comparison', 'Color', 'w');
	plot(k, mfSamp_awgn(k), '-', 'LineWidth', 1.5, 'Color', [0.10 0.45 0.75]); hold on;
	plot(k, mfSamp_h1(k), '-', 'LineWidth', 1.5, 'Color', [0.85 0.33 0.10]);
	plot(k, mfSamp_h2(k), '-', 'LineWidth', 1.5, 'Color', [0.20 0.65 0.20]);
	yline(0, 'k:');
	grid on;
	xlabel('Symbol Index');
	ylabel('MF Decision Sample');
	title(sprintf('MF Sample Stream Comparison @ Eb/N0 = %d dB', EYE_EBN0_DB));
	legend('AWGN', 'h1', 'h2', 'Location', 'best');
end

%% ---------------- Local Functions ----------------
function h = get_rrc_filter(sps, alpha, nTaps)
% Uses provided rrc_filter if available, otherwise designs an equivalent RRC.

if exist('rrc_filter', 'file') == 2
	h = rrc_filter(sps, alpha, nTaps);
	h = h(:).';
	return;
end

if exist('rcosdesign', 'file') ~= 2
	error(['Neither rrc_filter nor rcosdesign is available. ', ...
		   'Add the provided helper or Communications Toolbox.']);
end

% Convert desired taps to span symbols: length = span*sps + 1.
span = max(1, round((nTaps - 1) / sps));
h = rcosdesign(alpha, span, sps, 'sqrt');
end

function y = add_awgn_from_ebn0(x, EbN0dB, sps, k)
% Adds real AWGN using an Eb/N0 target for oversampled real baseband waveforms.

gamma_b = 10^(EbN0dB/10);

trim = min(10*sps, floor((numel(x)-1)/2));
if trim > 0
	xm = x(trim+1:end-trim);
else
	xm = x;
end

Px = mean(abs(xm).^2);
Es = Px * sps;
Eb = Es / k;
noiseVar = Eb / (2*gamma_b);

y = x + sqrt(noiseVar) * randn(size(x));
end

function plot_eye_custom(x, sps, nTraces, startIdx, rgb)
% Plots a 2-symbol eye by overlapping segments at symbol-rate offsets.

L = 2*sps;
N = numel(x);
t = (0:L-1)/sps;

maxTraces = floor((N - startIdx - L + 1) / sps);
nUse = max(0, min(nTraces, maxTraces));

hold on;
for n = 1:nUse
	idx0 = startIdx + (n-1)*sps;
	seg = x(idx0:idx0+L-1);
	plot(t, seg, '-', 'Color', rgb, 'LineWidth', 0.55);
end
end