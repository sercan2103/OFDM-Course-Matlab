clear all; close all; clc;

seed = randi([1, 10000]);
rng(seed);

N = 10000;
data = rand(1, N);

fprintf('=== Lab 1: Random Data Generation & Monte Carlo Simulation ===\n');
fprintf('Random Seed: %d (Change seed variable to get different results)\n\n', seed);
fprintf('Step 1: Generated %d uniform random samples\n', N);

fprintf('\nStep 2: Monte Carlo Simulation - 1 Million Coin Flips\n');
fprintf('----------------------------------------\n');

M = 1000000;
coin_flips = randi([0, 1], 1, M);

num_heads = sum(coin_flips);
num_tails = M - num_heads;

prob_heads_simulated = num_heads / M;
prob_tails_simulated = num_tails / M;

prob_theory = 0.5;

error_heads = abs(prob_heads_simulated - prob_theory) * 100;
error_tails = abs(prob_tails_simulated - prob_theory) * 100;

fprintf('Simulation Results (M = %d trials):\n', M);
fprintf('  Heads: %d (Probability = %.6f)\n', num_heads, prob_heads_simulated);
fprintf('  Tails: %d (Probability = %.6f)\n', num_tails, prob_tails_simulated);
fprintf('\nTheoretical Probability: %.6f\n', prob_theory);
fprintf('Error - Heads: %.4f%%\n', error_heads);
fprintf('Error - Tails: %.4f%%\n\n', error_tails);

figure('Position', [100 100 1200 800]);

subplot(2, 2, 1);
histogram(data, 50, 'Normalization', 'pdf', 'FaceColor', 'cyan', 'EdgeColor', 'black');
hold on;
yline(1, 'r-', 'LineWidth', 2, 'Label', 'Theoretical (p=1)');
title('Uniform Random Data Distribution (N=10,000)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Value'); 
ylabel('Probability Density');
legend('Histogram', 'Theory');
grid on;

subplot(2, 2, 2);
categories = {'Heads', 'Tails'};
counts = [num_heads, num_tails];
theoretical = [M*prob_theory, M*prob_theory];

x = 1:2;
bar(x, counts, 'FaceColor', 'cyan', 'EdgeColor', 'black', 'BarWidth', 0.6);
hold on;
plot(x, theoretical, 'ro-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Theoretical');

title('Monte Carlo: 1 Million Coin Flips', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Outcome');
ylabel('Count');
xticklabels(categories);
legend('Simulation', 'Theory');
grid on;
ylim([0 M*0.55]);

fprintf('Step 3: Convergence Analysis\n');
fprintf('----------------------------------------\n');

sample_sizes = [100, 1000, 10000, 100000, 1000000];
convergence_probs = [];

for sample_size = sample_sizes
    prob = sum(coin_flips(1:sample_size)) / sample_size;
    convergence_probs = [convergence_probs, prob];
    error = abs(prob - prob_theory) * 100;
    fprintf('N = %7d: P(Heads) = %.6f, Error = %.4f%%\n', sample_size, prob, error);
end

subplot(2, 2, 3);
semilogy(sample_sizes, abs(convergence_probs - prob_theory) * 100, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(sample_sizes, 1./sqrt(sample_sizes)*100, 'r--', 'LineWidth', 2, 'DisplayName', '1/√N (Theory)');
xlabel('Number of Trials');
ylabel('Error (%)', 'interpreter', 'tex');
title('Convergence Analysis: Error vs Sample Size', 'FontSize', 12, 'FontWeight', 'bold');
legend('Simulation Error', 'Theoretical Bound');
grid on;

subplot(2, 2, 4);
running_avg = cumsum(coin_flips) ./ (1:M);
sample_indices = round(logspace(2, 6, 100));
sample_indices = sample_indices(sample_indices <= M);
plot(sample_indices, running_avg(sample_indices), 'b-', 'LineWidth', 1.5);
hold on;
yline(prob_theory, 'r--', 'LineWidth', 2, 'DisplayName', 'Theoretical (0.5)');
set(gca, 'XScale', 'log');
xlabel('Number of Trials');
ylabel('Probability of Heads');
title('Running Average Convergence', 'FontSize', 12, 'FontWeight', 'bold');
ylim([0.45 0.55]);
legend('Running Average', 'Theory');
grid on;
