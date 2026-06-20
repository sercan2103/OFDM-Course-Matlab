%% Lab 3: Digital Modulation (BPSK, QPSK, 16-QAM)
% Implements tasks in lab3.md:
% 1) QPSK mod/demod
% 2) 16-QAM mod/demod
% 3) Constellation plots at SNR = 5, 10, 20 dB
% 4) BER comparison (BPSK, QPSK, 16-QAM) + theoretical curves
% 5) Sanity checks (normalization, Eb/N0 -> noise variance)
% Bonus) LLR examples

clear; clc; close all;
rng(42); % Reproducible results

fprintf('=== Lab 3: Digital Modulation Simulation ===\n');

%% Parameters
Nbits_demo = 2e5;      % For single-point demo (Task 1 & 2)
Nbits_ber  = 4e5;      % For BER sweep (Task 4)
EbN0_dB_demo = 10;     % Demo BER for Task 1 & 2
EbN0_dB_sweep = 0:16;  % Required range
scatterEbN0 = [5 10 20];

%% ------------------------ Task 1: QPSK ------------------------
% 1) Random bits
bits_qpsk = randi([0 1], Nbits_demo, 1);
if mod(numel(bits_qpsk),2) ~= 0
    bits_qpsk = [bits_qpsk; 0];
end

% 2) Serial-to-Parallel (pairs)
b2_qpsk = reshape(bits_qpsk, 2, []).';

% 3) I/Q mapping (signs from bits, then normalized to Es=1)
% bit 1 -> +1, bit 0 -> -1
I = 2*b2_qpsk(:,1) - 1;
Q = 2*b2_qpsk(:,2) - 1;
s_qpsk = (I + 1j*Q) / sqrt(2); % Es ~ 1

% 4) AWGN channel at Eb/N0
k_qpsk = 2;
EbN0_lin = 10^(EbN0_dB_demo/10);
N0_qpsk = 1 / (k_qpsk * EbN0_lin); % Es normalized to 1
sigma2_dim_qpsk = N0_qpsk/2;
n_qpsk = sqrt(sigma2_dim_qpsk) * (randn(size(s_qpsk)) + 1j*randn(size(s_qpsk)));
r_qpsk = s_qpsk + n_qpsk;

% 5) Threshold detector per axis
% (Equivalent decision boundaries on normalized constellation)
I_hat = real(r_qpsk) > 0;
Q_hat = imag(r_qpsk) > 0;

% 6) Parallel-to-Serial
bits_qpsk_hat = reshape([I_hat Q_hat].', [], 1);
bits_qpsk_ref = bits_qpsk(1:numel(bits_qpsk_hat));

% 7) BER
BER_qpsk_demo = mean(bits_qpsk_hat ~= bits_qpsk_ref);
fprintf('QPSK @ Eb/N0 = %2d dB -> BER(sim) = %.3e\n', EbN0_dB_demo, BER_qpsk_demo);

%% ------------------------ Task 2: 16-QAM ------------------------
% 1) Group bits in 4
bits_16qam = randi([0 1], Nbits_demo, 1);
rem4 = mod(numel(bits_16qam), 4);
if rem4 ~= 0
    bits_16qam = [bits_16qam; zeros(4-rem4,1)];
end
b4 = reshape(bits_16qam, 4, []).';

% 2) Map to 16-QAM grid levels {-3,-1,+1,+3} using Gray mapping on each axis
% Axis bits: [b1 b2] and [b3 b4]
Ilev = bits2pam4_gray(b4(:,1), b4(:,2));
Qlev = bits2pam4_gray(b4(:,3), b4(:,4));

% 3) Normalize by sqrt(10)
s_16qam = (Ilev + 1j*Qlev) / sqrt(10); % Es ~ 1

% 4) AWGN
k_16qam = 4;
N0_16qam = 1 / (k_16qam * EbN0_lin);
sigma2_dim_16qam = N0_16qam/2;
n_16qam = sqrt(sigma2_dim_16qam) * (randn(size(s_16qam)) + 1j*randn(size(s_16qam)));
r_16qam = s_16qam + n_16qam;

% 5) Demodulation / symbol decision
rI = real(r_16qam) * sqrt(10); % de-normalize to decision levels
rQ = imag(r_16qam) * sqrt(10);

[I_b1, I_b2] = pam4_gray2bits(rI);
[Q_b1, Q_b2] = pam4_gray2bits(rQ);

% 6) Symbols back to bits
bits_16qam_hat = reshape([I_b1 I_b2 Q_b1 Q_b2].', [], 1);
bits_16qam_ref = bits_16qam(1:numel(bits_16qam_hat));

% 7) BER
BER_16qam_demo = mean(bits_16qam_hat ~= bits_16qam_ref);
fprintf('16-QAM @ Eb/N0 = %2d dB -> BER(sim) = %.3e\n', EbN0_dB_demo, BER_16qam_demo);

%% ------------------ Task 3: Constellation Scatter Plots ------------------
Nsym_scatter = 4000;
figure('Name','Constellations at Different SNRs','Color','w','Position',[80 80 1200 700]);

for idx = 1:numel(scatterEbN0)
    snrDb = scatterEbN0(idx);
    snrLin = 10^(snrDb/10);

    % QPSK samples
    b = randi([0 1], 2*Nsym_scatter, 1);
    bp = reshape(b,2,[]).';
    s = ((2*bp(:,1)-1) + 1j*(2*bp(:,2)-1)) / sqrt(2);
    N0 = 1/(2*snrLin); % k=2
    n = sqrt(N0/2)*(randn(size(s))+1j*randn(size(s)));
    r = s + n;

    subplot(2,3,idx);
    plot(real(r), imag(r), '.', 'MarkerSize', 6);
    hold on;
    plot([-1 1]/sqrt(2), [-1 -1]/sqrt(2), 'ro', 'LineWidth', 1.2);
    plot([-1 1]/sqrt(2), [1 1]/sqrt(2), 'ro', 'LineWidth', 1.2);
    axis equal; grid on;
    xlabel('In-Phase'); ylabel('Quadrature');
    title(sprintf('QPSK, Eb/N0 = %d dB', snrDb));
    xlim([-2 2]); ylim([-2 2]);

    % 16-QAM samples
    b16 = randi([0 1], 4*Nsym_scatter, 1);
    g = reshape(b16,4,[]).';
    iL = bits2pam4_gray(g(:,1), g(:,2));
    qL = bits2pam4_gray(g(:,3), g(:,4));
    s16 = (iL + 1j*qL)/sqrt(10);
    N0_16 = 1/(4*snrLin); % k=4
    n16 = sqrt(N0_16/2)*(randn(size(s16))+1j*randn(size(s16)));
    r16 = s16 + n16;

    subplot(2,3,idx+3);
    plot(real(r16), imag(r16), '.', 'MarkerSize', 6);
    hold on;
    lv = [-3 -1 1 3]/sqrt(10);
    [xx,yy] = meshgrid(lv, lv);
    plot(xx(:), yy(:), 'ro', 'LineWidth', 1.2);
    axis equal; grid on;
    xlabel('In-Phase'); ylabel('Quadrature');
    title(sprintf('16-QAM, Eb/N0 = %d dB', snrDb));
    xlim([-1.5 1.5]); ylim([-1.5 1.5]);
end

%% --------------- Task 4: BER Curves (Sim + Theory) ----------------
BER_bpsk_sim  = zeros(size(EbN0_dB_sweep));
BER_qpsk_sim  = zeros(size(EbN0_dB_sweep));
BER_16q_sim   = zeros(size(EbN0_dB_sweep));
BER_bpsk_th   = zeros(size(EbN0_dB_sweep));
BER_qpsk_th   = zeros(size(EbN0_dB_sweep));
BER_16q_th    = zeros(size(EbN0_dB_sweep));

for ii = 1:numel(EbN0_dB_sweep)
    EbN0dB = EbN0_dB_sweep(ii);
    gamma_b = 10^(EbN0dB/10);

    %% BPSK
    b = randi([0 1], Nbits_ber, 1);
    s = 2*b - 1;
    sigma2 = 1/(2*gamma_b);       % per real dimension
    r = s + sqrt(sigma2)*randn(size(s));
    bh = r > 0;
    BER_bpsk_sim(ii) = mean(b ~= bh);
    BER_bpsk_th(ii) = qfunc(sqrt(2*gamma_b));

    %% QPSK (same BER as BPSK in AWGN with Gray mapping)
    bq = randi([0 1], Nbits_ber, 1);
    if mod(numel(bq),2)~=0, bq = [bq;0]; end
    p = reshape(bq,2,[]).';
    sq = ((2*p(:,1)-1) + 1j*(2*p(:,2)-1))/sqrt(2);
    N0 = 1/(2*gamma_b);           % k=2
    nq = sqrt(N0/2)*(randn(size(sq))+1j*randn(size(sq)));
    rq = sq + nq;
    bI = real(rq)>0;
    bQ = imag(rq)>0;
    bhq = reshape([bI bQ].', [], 1);
    BER_qpsk_sim(ii) = mean(bhq ~= bq(1:numel(bhq)));
    BER_qpsk_th(ii)  = qfunc(sqrt(2*gamma_b));

    %% 16-QAM (Gray)
    b16 = randi([0 1], Nbits_ber, 1);
    m = mod(numel(b16),4);
    if m~=0, b16 = [b16; zeros(4-m,1)]; end
    g = reshape(b16,4,[]).';
    I = bits2pam4_gray(g(:,1), g(:,2));
    Q = bits2pam4_gray(g(:,3), g(:,4));
    s16 = (I + 1j*Q)/sqrt(10);
    N0_16 = 1/(4*gamma_b);        % k=4
    n16 = sqrt(N0_16/2)*(randn(size(s16))+1j*randn(size(s16)));
    r16 = s16 + n16;

    rI = real(r16)*sqrt(10);
    rQ = imag(r16)*sqrt(10);
    [i1,i2] = pam4_gray2bits(rI);
    [q1,q2] = pam4_gray2bits(rQ);
    bh16 = reshape([i1 i2 q1 q2].', [], 1);
    BER_16q_sim(ii) = mean(bh16 ~= b16(1:numel(bh16)));

    % Common analytical approximation for Gray-coded square 16-QAM
    BER_16q_th(ii) = (4/4)*(1-1/sqrt(16))*qfunc(sqrt((3*4/(16-1))*gamma_b));

    fprintf('Eb/N0=%2d dB | BER BPSK=%.3e, QPSK=%.3e, 16QAM=%.3e\n', ...
        EbN0dB, BER_bpsk_sim(ii), BER_qpsk_sim(ii), BER_16q_sim(ii));
end

figure('Name','BER Comparison','Color','w');
semilogy(EbN0_dB_sweep, BER_bpsk_sim, 'bo-', 'LineWidth',1.4, 'MarkerSize',5); hold on;
semilogy(EbN0_dB_sweep, BER_qpsk_sim, 'rs-', 'LineWidth',1.4, 'MarkerSize',5);
semilogy(EbN0_dB_sweep, BER_16q_sim, 'md-', 'LineWidth',1.4, 'MarkerSize',5);
semilogy(EbN0_dB_sweep, BER_bpsk_th, 'b--', 'LineWidth',1.8);
semilogy(EbN0_dB_sweep, BER_qpsk_th, 'r--', 'LineWidth',1.8);
semilogy(EbN0_dB_sweep, BER_16q_th, 'm--', 'LineWidth',1.8);
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BER Comparison: BPSK vs QPSK vs 16-QAM (AWGN)');
legend('BPSK Sim','QPSK Sim','16-QAM Sim','BPSK Theory','QPSK Theory','16-QAM Theory', ...
    'Location','southwest');
ylim([1e-5 1]);

%% ------------------------ Task 5: Sanity Checks ------------------------
Es_qpsk = mean(abs(s_qpsk).^2);
Es_16   = mean(abs(s_16qam).^2);

fprintf('\n=== Sanity Checks ===\n');
fprintf('Avg symbol energy QPSK:   Es = %.4f (target ~1)\n', Es_qpsk);
fprintf('Avg symbol energy 16-QAM: Es = %.4f (target ~1)\n', Es_16);
fprintf('QPSK: sigma^2 per dim from Eb/N0=%d dB -> %.4e\n', EbN0_dB_demo, sigma2_dim_qpsk);
fprintf('16QAM: sigma^2 per dim from Eb/N0=%d dB -> %.4e\n', EbN0_dB_demo, sigma2_dim_16qam);

%% ---------------------- Bonus: LLR Computation ----------------------
% A) Manual BPSK LLR: LLR = 2r/sigma^2
Nb_llr = 6e4;
b = randi([0 1], Nb_llr, 1);
s = 2*b - 1;
gamma_b = 10^(6/10);      % Example Eb/N0 = 6 dB
sigma2 = 1/(2*gamma_b);
r = s + sqrt(sigma2)*randn(size(s));
llr_bpsk = 2*r/sigma2;

figure('Name','Bonus - LLR Histograms','Color','w');
subplot(1,2,1);
histogram(llr_bpsk(b==0), 80, 'Normalization','pdf'); hold on;
histogram(llr_bpsk(b==1), 80, 'Normalization','pdf');
grid on;
xlabel('LLR'); ylabel('PDF');
title('BPSK LLR Distribution (Eb/N0=6 dB)');
legend('bit=0','bit=1');

% B) 16-QAM soft-demod LLR (if toolbox available)
subplot(1,2,2);
if exist('qamdemod','file') == 2
    M = 16;
    Nsym = 3e4;
    bits = randi([0 1], Nsym*log2(M), 1);

    % Use built-in to generate Gray-coded normalized symbols
    s_idx = bi2de(reshape(bits, log2(M), []).', 'left-msb');
    s = qammod(s_idx, M, 'gray', 'UnitAveragePower', true);

    EbN0_dB = 10;
    gamma_b = 10^(EbN0_dB/10);
    N0 = 1/(log2(M)*gamma_b);
    n = sqrt(N0/2)*(randn(size(s)) + 1j*randn(size(s)));
    r = s + n;

    llr = qamdemod(r, M, 'gray', 'UnitAveragePower', true, ...
        'OutputType', 'llr', 'NoiseVariance', N0);

    histogram(llr(:), 100, 'Normalization','pdf');
    title('16-QAM LLR Histogram (Toolbox)');
    xlabel('LLR'); ylabel('PDF'); grid on;
else
    axis off;
    text(0.05, 0.5, 'qamdemod not found (Communications Toolbox missing).', 'FontSize', 11);
    title('16-QAM LLR Histogram');
end

fprintf('\nDone. Generated: constellation plots, BER curve, sanity checks, and LLR figures.\n');

%% ---------------------------- Local Functions ----------------------------
function a = bits2pam4_gray(b1, b2)
% Gray mapping: 00->-3, 01->-1, 11->+1, 10->+3
idx = 2*b1 + b2; % 0,1,2,3 for 00,01,10,11
a = zeros(size(idx));
a(idx==0) = -3; % 00
a(idx==1) = -1; % 01
a(idx==3) = +1; % 11
a(idx==2) = +3; % 10
end

function [b1, b2] = pam4_gray2bits(x)
% Decision thresholds at -2, 0, +2 for levels -3,-1,+1,+3
xq = zeros(size(x));
xq(x < -2) = -3;
xq(x >= -2 & x < 0) = -1;
xq(x >= 0 & x < 2) = +1;
xq(x >= 2) = +3;

b1 = zeros(size(xq));
b2 = zeros(size(xq));

% Inverse Gray map:
% -3 -> 00
% -1 -> 01
% +1 -> 11
% +3 -> 10
b1(xq==+1 | xq==+3) = 1;
b2(xq==-1 | xq==+1) = 1;
end
