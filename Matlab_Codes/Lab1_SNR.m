%% Lab 1: AWGN & SNR
% TODO-1: Generate clean signal
fs=8000; t=0:1/fs:1; f0=440;
s = sin(2*pi*f0*t); % Fill in!
% TODO-2: Signal power
P_s = mean(s.^2);
% TODO-3: Sigma from SNR_dB=10
SNR_lin = 10^(10/10);
sigma = sqrt(P_s/SNR_lin);
% TODO-4: Add AWGN
n = sigma*randn(size(t)); r = s + n;
% TODO-5: Verify SNR
P_n = mean(n.^2);
SNR_chk = 10*log10(P_s/P_n);
fprintf('%.1f dB\n', SNR_chk);
% TODO-6: Plot 6-panel subplot
figure('Position', [100 100 1200 600]);

% Panel 1: Clean Signal
subplot(2,3,1);
plot(t(1:160), s(1:160), 'b', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Amplitude');
title('Clean Signal');
grid on;

% Panel 2: Noise Histogram
subplot(2,3,2);
n_noise = sigma*randn(1, 8000);
histogram(n_noise, 50, 'FaceColor', 'r', 'EdgeColor', 'k');
xlabel('Amplitude');
ylabel('Frequency');
title('Noise Histogram');
grid on;

% Panel 3: SNR = +20 dB
snr_vals = [20, 10, 0, -5];
titles = {'VERY GOOD (+20 dB)', 'GOOD (+10 dB)', 'BOUNDARY (0 dB)', 'POOR (-5 dB)'};
for i = 1:4
    SNR_dB = snr_vals(i);
    SNR_lin = 10^(SNR_dB/10);
    sigma_i = sqrt(P_s/SNR_lin);
    n_i = sigma_i*randn(size(t));
    r_i = s + n_i;
    
    subplot(2,3,i+2);
    plot(t(1:160), r_i(1:160), 'LineWidth', 1);
    xlabel('Time (ms)');
    ylabel('Amplitude');
    title(sprintf('SNR = %+d dB', SNR_dB));
    grid on;
end
% TODO-7: sound() at each SNR
for i = 1:length(snr_vals)
    SNR_dB = snr_vals(i);
    SNR_lin = 10^(SNR_dB/10);
    sigma_i = sqrt(P_s/SNR_lin);
    n_i = sigma_i*randn(size(t));
    r_i = s + n_i;
    
    fprintf('Playing SNR = %d dB. Press any key to continue...\n', SNR_dB);
    sound(r_i, fs);
    pause;
end


%Bonus 1: Plot SNR accuracy vs N samples (100 to 100,000)
sample_sizes = [100, 500, 1000, 5000, 10000, 50000, 100000];
snr_errors = zeros(size(sample_sizes));
for i = 1:length(sample_sizes)
    N = sample_sizes(i);
    s_sample = sin(2*pi*f0*(0:N-1)/fs);
    P_s_sample = mean(s_sample.^2);
    
    SNR_dB_target = 10; % Target SNR in dB
    SNR_lin_target = 10^(SNR_dB_target/10);
    sigma_sample = sqrt(P_s_sample/SNR_lin_target);
    
    n_sample = sigma_sample*randn(size(s_sample));
    r_sample = s_sample + n_sample;
    
    P_n_sample = mean(n_sample.^2);
    SNR_dB_estimated = 10*log10(P_s_sample/P_n_sample);
    
    snr_errors(i) = abs(SNR_dB_estimated - SNR_dB_target);
end
figure;
semilogx(sample_sizes, snr_errors, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of Samples (N)');            
ylabel('SNR Estimation Error (dB)');
title('SNR Estimation Error vs Number of Samples', 'FontSize', 12, 'FontWeight', 'bold');
grid on;    

%Bonus 2:
%Colored noise — compare PSD with pwelch()

figure;
fs = 8000;
t = 0:1/fs:1;
f0 = 440;
s = sin(2*pi*f0*t);

% Generate colored noise (e.g., pink noise)
n_colored = filter([1; -0.9], 1, randn(size(t))'); % Simple pink noise generation
n_colored = n_colored / max(abs(n_colored)); % Normalize

% Add colored noise to signal
r_colored = s + 0.5 * n_colored; % Scale colored noise

% Compare PSD using pwelch
[pxx, f] = pwelch(r_colored, 512, [], [], fs);
[pxx_s, ~] = pwelch(s, 512, [], [], fs);
[pxx_n, ~] = pwelch(n_colored, 512, [], [], fs);

figure('Position', [100 100 1000 600]);
semilogy(f, pxx, 'b-', 'LineWidth', 2);
hold on;
semilogy(f, pxx_s, 'r--', 'LineWidth', 2);
semilogy(f, pxx_n, 'g:', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density');
title('PSD Comparison: Colored Noise vs Signal vs Noise');
legend('Colored Noise + Signal (blue solid)', 'Signal Only (red dashed)', 'Colored Noise Only (green dotted)', 'Location', 'best');
grid on;
hold off;

%Bonus 3:
%Estimate SNR WITHOUT knowing original signal
%Assume we only have the noisy signal r_colored
P_r = mean(r_colored.^2);
P_n_est = mean((r_colored - mean(r_colored)).^2); % Estimate noise power
SNR_est = 10*log10(P_r/P_n_est);
fprintf('Estimated SNR = %.1f dB\n', SNR_est);
