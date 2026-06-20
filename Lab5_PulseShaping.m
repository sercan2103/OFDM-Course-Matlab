%% Lab 5: Pulse Shaping and Matched Filtering
% 1) Generate BPSK symbols
% 2) Upsample and RC/RRC pulse shaping
% 3) Add AWGN noise
% 4) Apply matched filter at Rx
% 5) Plot eye diagrams
% 6) BER with MF vs without MF
% 7) Show matched filter improves BER

clear; clc; close all;
rng(45);

fprintf('=== Lab 5: Pulse Shaping and Matched Filtering ===\n');

%% Parameters
Nsym       = 2e5;         % Number of symbols for BER sweep
sps        = 8;           % Samples per symbol
span       = 10;          % RRC span in symbols
alpha      = 0.35;        % Roll-off factor (main case)
EbN0_dB    = 0:2:10;      % BER sweep points

% For eye diagram experiments
eyeNsym    = 4000;
eyeAlpha   = [0.2 0.35 0.8];
eyeSNRdB   = [12 8 4];

%% Filter design (requires Communications Toolbox)
if exist('rcosdesign','file') ~= 2
	error(['Function rcosdesign not found. Please enable/install Communications Toolbox ', ...
		   'to run Lab 5.']);
end

h_rrc = rcosdesign(alpha, span, sps, 'sqrt');   % Tx/Rx matched pair
h_rc  = rcosdesign(alpha, span, sps, 'normal'); % For spectrum comparison
h_rect = ones(1, sps);                           % Rectangular pulse

gd = (length(h_rrc)-1)/2; % Group delay of one RRC filter in samples

%% 1) Generate BPSK symbols
bits = randi([0 1], Nsym, 1);
sym  = 2*bits - 1; % 0->-1, 1->+1

%% 2) Upsample + pulse shaping (RRC)
tx_rrc = upfirdn(sym, h_rrc, sps, 1);

%% 3) Add AWGN + 4) Matched filter + 5) Eye diagrams
figure('Name','Eye Diagrams (alpha and SNR sweeps)','Color','w','Position',[100 100 1200 700]);

% Top row: alpha sweep, fixed SNR
fixedSNR = 8;
for k = 1:numel(eyeAlpha)
	a = eyeAlpha(k);
	htx = rcosdesign(a, span, sps, 'sqrt');
	hrx = htx;
	gd_k = (length(htx)-1)/2;

	b = randi([0 1], eyeNsym, 1);
	s = 2*b - 1;
	tx = upfirdn(s, htx, sps, 1);

	rx = addAwgnFromEbN0(tx, fixedSNR, sps, 1);
	y  = conv(rx, hrx, 'full');

	subplot(2,3,k);
	plotEyeCustom(y, sps, 250, 2*gd_k + 1, [0.30 0.60 0.80], 0.08);
	grid on; xlabel('Time (T)'); ylabel('Amplitude');
	title(sprintf('Eye: \\alpha=%.2f, Eb/N0=%d dB', a, fixedSNR));
	xlim([0 2]);
end

% Bottom row: SNR sweep, fixed alpha
fixedAlpha = alpha;
htx = rcosdesign(fixedAlpha, span, sps, 'sqrt');
hrx = htx;
gd_eye = (length(htx)-1)/2;

for k = 1:numel(eyeSNRdB)
	snrDb = eyeSNRdB(k);

	b = randi([0 1], eyeNsym, 1);
	s = 2*b - 1;
	tx = upfirdn(s, htx, sps, 1);

	rx = addAwgnFromEbN0(tx, snrDb, sps, 1);
	y  = conv(rx, hrx, 'full');

	subplot(2,3,k+3);
	plotEyeCustom(y, sps, 250, 2*gd_eye + 1, [0.10 0.45 0.75], 0.08);
	grid on; xlabel('Time (T)'); ylabel('Amplitude');
	title(sprintf('Eye: \\alpha=%.2f, Eb/N0=%d dB', fixedAlpha, snrDb));
	xlim([0 2]);
end

sgtitle('Eye diagrams: higher \\alpha and lower Eb/N0 close the eye');

%% Spectrum comparison: Rectangular vs RC(alpha=0.35)
Nspec = 2e4;
b_spec = randi([0 1], Nspec, 1);
s_spec = 2*b_spec - 1;

tx_rect = upfirdn(s_spec, h_rect, sps, 1);
tx_rc   = upfirdn(s_spec, h_rc,   sps, 1);

[f_rect, P_rect_dB] = estimateSpectrumDb(tx_rect, sps);
[f_rc,   P_rc_dB]   = estimateSpectrumDb(tx_rc,   sps);

figure('Name','Spectrum Comparison','Color','w');
plot(f_rect, P_rect_dB, 'Color', [0.90 0.45 0.40], 'LineWidth', 0.9); hold on;
plot(f_rc,   P_rc_dB,   'Color', [0.10 0.45 0.75], 'LineWidth', 1.7);
grid on;
xlabel('Frequency (1/T)');
ylabel('Power (dB)');
title('Spectrum Comparison: Rectangular vs RC(\\alpha=0.35)');
legend('Rect', 'RC \\alpha=0.35', 'Location', 'southwest');
xlim([-2 2]);
ylim([-70 5]);

%% 6) BER with matched filter vs without matched filter
ber_no_mf = zeros(size(EbN0_dB));
ber_mf    = zeros(size(EbN0_dB));
ber_th    = zeros(size(EbN0_dB));

for ii = 1:numel(EbN0_dB)
	ebn0 = EbN0_dB(ii);

	b = randi([0 1], Nsym, 1);
	s = 2*b - 1;
	tx = upfirdn(s, h_rrc, sps, 1);

	% 3) AWGN at waveform level
	rx = addAwgnFromEbN0(tx, ebn0, sps, 1);

	% Without matched filter (sample directly after Tx pulse-shaping + noise)
	start_no_mf = gd + 1;
	z_no_mf = rx(start_no_mf:sps:start_no_mf + (Nsym-1)*sps);
	bh_no_mf = z_no_mf > 0;
	ber_no_mf(ii) = mean(bh_no_mf ~= b);

	% 4) With matched filter
	y_mf = conv(rx, h_rrc, 'full');
	start_mf = 2*gd + 1;
	z_mf = y_mf(start_mf:sps:start_mf + (Nsym-1)*sps);
	bh_mf = z_mf > 0;
	ber_mf(ii) = mean(bh_mf ~= b);

	gamma_b = 10^(ebn0/10);
	ber_th(ii) = qfunc(sqrt(2*gamma_b));

	fprintf('Eb/N0=%2d dB | BER(no MF)=%.3e | BER(MF)=%.3e | Theory=%.3e\n', ...
		ebn0, ber_no_mf(ii), ber_mf(ii), ber_th(ii));
end

%% 7) "Proof": matched filter improves BER
figure('Name','BER: MF vs No-MF','Color','w');
semilogy(EbN0_dB, ber_th,    'k--', 'LineWidth', 1.8); hold on;
semilogy(EbN0_dB, ber_mf,    'o-',  'Color', [0.10 0.45 0.75], 'LineWidth', 2.0, 'MarkerSize', 7);
semilogy(EbN0_dB, ber_no_mf, 's-',  'Color', [0.85 0.33 0.10], 'LineWidth', 2.0, 'MarkerSize', 7);
grid on;
xlabel('Eb/N0 (dB)');
ylabel('BER');
title('BER vs Eb/N0: Matched Filter Gain');
legend('Theory (BPSK)', 'Simulation (with MF)', 'Simulation (without MF)', ...
	'Location', 'southwest');
ylim([1e-5 1]);

improvedAtAll = all(ber_mf <= ber_no_mf);
avgGainRatio = mean(ber_no_mf ./ max(ber_mf, eps));

fprintf('\n=== Matched Filter Conclusion ===\n');
fprintf('MF better or equal at all tested Eb/N0 points: %d\n', improvedAtAll);
fprintf('Average BER improvement factor (No-MF / MF): %.2f x\n', avgGainRatio);
fprintf('Done.\n');

%% ---------------- Local Functions ----------------
function y = addAwgnFromEbN0(x, EbN0dB, sps, k)
% Adds real AWGN to a real waveform x using Eb/N0 mapping.
% k = bits per symbol (BPSK -> k=1)

gamma_b = 10^(EbN0dB/10);

% Estimate Es from waveform power (trim filter transients for better estimate)
trim = min(10*sps, floor((numel(x)-1)/2));
if trim > 0
	xm = x(trim+1:end-trim);
else
	xm = x;
end

Px = mean(abs(xm).^2);   % Average power at sample rate Fs = sps/T
Es = Px * sps;           % Symbol energy estimate (for unit-energy pulse shaping)
Eb = Es / k;

% Real baseband AWGN variance per sample
noiseVar = Eb / (2*gamma_b);
y = x + sqrt(noiseVar)*randn(size(x));
end

function plotEyeCustom(x, sps, nTraces, startIdx, rgb, alphaVal)
% Plots a 2-symbol eye diagram using overlapped traces.

L = 2*sps;
N = length(x);
t = (0:L-1)/sps;
hold on;

maxTraces = floor((N - startIdx - L + 1)/sps);
nUse = min(nTraces, maxTraces);

for n = 1:nUse
	idx0 = startIdx + (n-1)*sps;
	seg = x(idx0:idx0+L-1);
	plot(t, seg, '-', 'Color', [rgb alphaVal], 'LineWidth', 0.8);
end
end

function [f, PdB] = estimateSpectrumDb(x, Fs)
% Simple FFT-based spectrum estimate normalized to 0 dB peak.

Nfft = 2^nextpow2(min(length(x), 65536));
X = fftshift(fft(x, Nfft));
P = abs(X).^2;
P = P ./ max(P + eps);

f = (-Nfft/2:Nfft/2-1)*(Fs/Nfft);
PdB = 10*log10(P + eps);

% light smoothing for cleaner visual comparison
win = 9;
PdB = conv(PdB, ones(1,win)/win, 'same');
end