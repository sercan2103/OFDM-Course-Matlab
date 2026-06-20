% === Op-count estimator via calibration ===

% Step 1: calibrate MACs/sec on this machine
N = 1e6;
a = randn(N, 1);
b = randn(N, 1);
c = randn(N, 1);
% warm up JIT
for i = 1:5, d = a.*b + c; end

% time 10 runs of N multiply-adds
t0 = tic;
for i = 1:10
    d = a.*b + c;   % N multiply-accumulates
end
t_per_mac = toc(t0) / (10 * N);   % seconds per MAC
macs_per_sec = 1 / t_per_mac;
fprintf('Calibration: %.2e MACs/sec\n', macs_per_sec);

% Step 2: algorithmic frame time (subtract rng overhead from profiler)
t_frame_total    = 0.091;   % profiler total
t_rng_overhead   = 0.035;   % rng inside local_pilot_values (not algorithm)
t_frame_algo     = t_frame_total - t_rng_overhead;  % 0.056 s

% Step 3: estimated op-count
est_ops = round(t_frame_algo * macs_per_sec);
fprintf('Estimated op-count: %.2e  (~%d MACs/frame)\n', est_ops, est_ops);