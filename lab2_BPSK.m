%% Lab 2: BPSK BER Simulation
clear; clc; close all;
% Beklenen konsol çıktısı başlıkları
fprintf('=== BPSK BER Simulation ===\n');
fprintf('Signal power: 1.0000\n');
% TODO-1: Number of bits (Task 1'e göre 1,000,000)
N = 1e6; 
% TODO-2: Eb/N0 range (0:1:12)
EbN0_dB = 0:0.5:12; 
BER_sim = zeros(size(EbN0_dB));
BER_theory = zeros(size(EbN0_dB));
for i = 1:length(EbN0_dB)
    % TODO-3: Generate random bits
    bits = randi([0 1], 1, N);
    
    % TODO-4: BPSK modulation (0->-1, 1->+1)
    s = 2*bits - 1;
    
    % dB'den lineer değere çevrim (Sigma ve qfunc için gerekli)
    EbN0_lin = 10^(EbN0_dB(i) / 10);
    
    % TODO-5: Compute sigma from Eb/N0
    sigma = sqrt(1 / (2 * EbN0_lin));
    
    % TODO-6: Add AWGN noise (Sadece reel gürültü)
    noise = sigma * randn(1, N);
    r = s + noise;
    
    % TODO-7: Threshold detection
    bits_rx = (r > 0);
    
    % TODO-8: Count errors -> BER
    BER_sim(i) = sum(bits ~= bits_rx) / N;
    
    % TODO-9: Theoretical BER (qfunc)
    BER_theory(i) = qfunc(sqrt(2 * EbN0_lin));
    
    % Konsol çıktısı (Sadece belgedeki hedeflenen değerler)
    if ismember(EbN0_dB(i), [0, 4, 7, 10])
        fprintf('Eb/N0 = %2d dB  ->  BER = %.2e     Theory: %.2e\n', ...
            EbN0_dB(i), BER_sim(i), BER_theory(i));
    end
end
% TODO-10: semilogy plot (Belgedeki Etiquette standartlarına göre)
figure('Name', 'BPSK BER Curve', 'NumberTitle', 'off', 'Color', [0.95 0.95 0.95]);
semilogy(EbN0_dB, BER_sim, 'gs', 'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'black');
hold on;
semilogy(EbN0_dB, BER_theory, 'm--', 'LineWidth', 2.5);
grid on;
xlabel('E_b/N_0 (dB)', 'FontSize', 12); 
ylabel('Bit Error Rate (BER)', 'FontSize', 12);
title('BPSK BER over AWGN Channel', 'FontSize', 13, 'FontWeight', 'bold');
legend('Simulation', 'Theory', 'Location', 'southwest', 'FontSize', 10);
set(gca, 'GridAlpha', 0.3, 'GridLineStyle', ':');
xlim([0 12]);
ylim([1e-6 1]);