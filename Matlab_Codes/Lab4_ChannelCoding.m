%% Lab 4: Progressive Channel Coding BER (QPSK/BPSK-equivalent)
% 1) Uncoded baseline (same BER as QPSK Gray in AWGN)
% 2) Repetition-3 (encode + majority vote)
% 3) Convolutional code R=1/2, K=7 (convenc)
% 4) Hard-decision Viterbi BER curve
% 5) Soft-decision Viterbi BER curve
% 6) Plot all 5 curves on one figure
% Bonus) Punctured R=2/3 curve

clear; clc; close all;
rng(44);

fprintf('=== Lab 4: Progressive Coding (Rep -> Conv Hard -> Conv Soft) ===\n');

%% Parameters
Nbits       = 2e5;          % Information bits per Eb/N0 point
EbN0_dB     = 0:0.5:12;     % Sweep
targetBER   = 1e-5;         % For coding gain table

% Convolutional code (R=1/2, K=7)
K       = 7;
trellis = poly2trellis(K, [171 133]);
tblen   = 5*(K-1);

% Tool check for conv/Viterbi tasks
hasConvToolbox = (exist('convenc','file') == 2) && (exist('vitdec','file') == 2) && (exist('poly2trellis','file') == 2);
if ~hasConvToolbox
    error(['Communications Toolbox functions not found (convenc/vitdec/poly2trellis). ', ...
           'Install/enable toolbox to run Lab 4 convolutional parts.']);
end

%% BER arrays
ber_uncoded   = zeros(size(EbN0_dB));
ber_rep3      = zeros(size(EbN0_dB));
ber_conv_hard = zeros(size(EbN0_dB));
ber_conv_soft = zeros(size(EbN0_dB));
ber_theory    = zeros(size(EbN0_dB));
ber_punc23    = nan(size(EbN0_dB)); % Bonus

puncpat23 = [1 1 0 1]; % Mother R=1/2 -> punctured R=2/3

%% Main sweep
for ii = 1:numel(EbN0_dB)
    ebn0dB = EbN0_dB(ii);
    gamma_b = 10^(ebn0dB/10);

    % Shared source bits per point
    b = randi([0 1], Nbits, 1);

    %% 1) Uncoded baseline (BPSK BER = QPSK Gray BER)
    sigma = sqrt(1/(2*gamma_b));
    x = 1 - 2*b;                          % 0->+1, 1->-1
    r = x + sigma*randn(size(x));
    bh = r < 0;
    ber_uncoded(ii) = mean(b ~= bh);

    %% 2) Repetition-3 + majority vote (rewritten)
    % Encode: [b1 b2 ...] -> [b1 b1 b1 b2 b2 b2 ...]
    b_rep = zeros(3*Nbits, 1, 'like', b);
    b_rep(1:3:end) = b;
    b_rep(2:3:end) = b;
    b_rep(3:3:end) = b;

    sigma = sqrt(1/(2*gamma_b));
    x_rep = 1 - 2*b_rep;                  % BPSK: 0->+1, 1->-1
    r_rep = x_rep + sigma*randn(size(x_rep));

    % Hard decision on coded bits
    bh_rep_coded = double(r_rep < 0);

    % Decode by majority vote per triplet
    bh_rep_triplets = reshape(bh_rep_coded, 3, Nbits);
    bh_rep = (sum(bh_rep_triplets, 1) >= 2).';

    ber_rep3(ii) = mean(bh_rep ~= b);

    %% 3+4) Convolutional R=1/2 + hard-decision Viterbi
    c = convenc(b, trellis);
    sigma = sqrt(1/(2*gamma_b));
    x_c = 1 - 2*c;
    r_c = x_c + sigma*randn(size(x_c));

    c_hard = r_c < 0;
    b_hat_hard = vitdec(c_hard, trellis, tblen, 'trunc', 'hard');
    ber_conv_hard(ii) = mean(b_hat_hard ~= b);

    %% 5) Soft-decision Viterbi (unquantized)
    % For 'unquant', positive metric favors bit 0 and negative favors bit 1
    b_hat_soft = vitdec(r_c, trellis, tblen, 'trunc', 'unquant');
    ber_conv_soft(ii) = mean(b_hat_soft ~= b);

    %% Theory (BPSK/QPSK-Gray)
    ber_theory(ii) = qfunc(sqrt(2*gamma_b));

    %% Bonus) Punctured R=2/3 (hard-decision)
    c23 = convenc(b, trellis, puncpat23);
    sigma = sqrt(1/(2*gamma_b));
    x23 = 1 - 2*c23;
    r23 = x23 + sigma*randn(size(x23));
    c23_hard = r23 < 0;
    b23_hat = vitdec(c23_hard, trellis, tblen, 'trunc', 'hard', puncpat23);
    ber_punc23(ii) = mean(b23_hat ~= b);

    fprintf('Eb/N0=%4.1f dB | Uncoded=%.3e | Rep3=%.3e | ConvH=%.3e | ConvS=%.3e\n', ...
        ebn0dB, ber_uncoded(ii), ber_rep3(ii), ber_conv_hard(ii), ber_conv_soft(ii));
end

%% 6) Plot 5 required BER curves on one graph
figure('Name','Lab 4 Expected BER Output','Color','w');
semilogy(EbN0_dB, ber_uncoded,   '-',  'Color',[0.55 0.55 0.55], 'LineWidth',2.0); hold on;
semilogy(EbN0_dB, ber_rep3,      '-',  'Color',[0.85 0.33 0.30], 'LineWidth',2.2);
semilogy(EbN0_dB, ber_conv_hard, '-',  'Color',[0.20 0.55 0.85], 'LineWidth',2.2);
semilogy(EbN0_dB, ber_conv_soft, '-',  'Color',[0.00 0.60 0.60], 'LineWidth',2.2);
semilogy(EbN0_dB, ber_theory,    '--', 'Color',[0.70 0.70 0.70], 'LineWidth',1.8);

if all(~isnan(ber_punc23))
    semilogy(EbN0_dB, ber_punc23, '-', 'Color',[0.20 0.75 0.35], 'LineWidth',2.0);
end

grid on;
xlabel('Eb/N0 (dB)');
ylabel('BER');
title('Lab 4 Expected BER Output');
legendItems = {'QPSK Uncoded','Rep-3 (R=1/3)','Conv Hard','Conv Soft','Theory (BPSK)'};
if all(~isnan(ber_punc23))
    legendItems{end+1} = 'Punctured Conv (R=2/3)';
end
legend(legendItems, 'Location','southwest');
ylim([1e-7 1e0]);
xlim([min(EbN0_dB) max(EbN0_dB)]);

%% Deliverable 2: coding gain table @ BER = 1e-5
EbN0_u = findEbN0AtTarget(EbN0_dB, ber_uncoded, targetBER);
EbN0_r = findEbN0AtTarget(EbN0_dB, ber_rep3, targetBER);
EbN0_h = findEbN0AtTarget(EbN0_dB, ber_conv_hard, targetBER);
EbN0_s = findEbN0AtTarget(EbN0_dB, ber_conv_soft, targetBER);

fprintf('\n=== Coding Gain Table at BER = %.1e ===\n', targetBER);
fprintf('Uncoded baseline : Eb/N0 = %6.3f dB\n', EbN0_u);
fprintf('Rep-3 decoder    : Eb/N0 = %6.3f dB | Gain = %6.3f dB\n', EbN0_r, EbN0_u - EbN0_r);
fprintf('Conv Hard decoder: Eb/N0 = %6.3f dB | Gain = %6.3f dB\n', EbN0_h, EbN0_u - EbN0_h);
fprintf('Conv Soft decoder: Eb/N0 = %6.3f dB | Gain = %6.3f dB\n', EbN0_s, EbN0_u - EbN0_s);

if all(~isnan(ber_punc23))
    EbN0_p = findEbN0AtTarget(EbN0_dB, ber_punc23, targetBER);
    fprintf('Punctured R=2/3  : Eb/N0 = %6.3f dB | Gain = %6.3f dB\n', EbN0_p, EbN0_u - EbN0_p);
end

fprintf('\nDone. Curves should show progressive gain: rep -> conv hard -> conv soft.\n');

%% Local function
function xAtTarget = findEbN0AtTarget(xdB, ber, target)
% Returns Eb/N0(dB) where BER reaches target (log-domain interpolation)
% If target is not bracketed by simulated BER, returns NaN.
ber = ber(:); xdB = xdB(:);

if any(ber <= target) && any(ber >= target)
    i2 = find(ber <= target, 1, 'first');
    if i2 == 1
        xAtTarget = xdB(1);
        return;
    end
    i1 = i2 - 1;

    y1 = log10(ber(i1));
    y2 = log10(ber(i2));
    yt = log10(target);

    if abs(y2 - y1) < eps
        xAtTarget = xdB(i2);
    else
        xAtTarget = xdB(i1) + (xdB(i2)-xdB(i1)) * (yt-y1)/(y2-y1);
    end
else
    xAtTarget = NaN;
end
end
