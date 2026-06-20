% VERIFY_DOPPLER0  Impairment-isolation diagnostic.
%
% Goal: figure out which impairment is responsible for the ~0.1-0.4 BER we
% see with the pass-through stubs. We toggle Doppler / STO / CFO off one
% at a time (and all together) at a single Eb/N0 and report the BER.
%
% Expected outcome (with the stubs unchanged):
%   * ALL impairments off  -> BER at the AWGN floor (proves the chain is sound)
%   * only Doppler off     -> still ~0.3-0.4 (CFO/STO dominate)
%   * only STO off         -> still bad (CFO dominates)
%   * only CFO off         -> noticeably better (now only STO + Doppler hurt)
%   * baseline (all on)    -> ~0.3-0.4 (random output)
%
% If the "ALL OFF" row is NOT at the AWGN floor, the chain itself is buggy
% -- look at channel_estimate / equalize / qpsk_demap.
% If "ALL OFF" IS at the floor, the chain is fine and the only fix is to
% implement sto_estimate_correct() and cfo_estimate_correct().
%
% Note: with the CP-self-correlation baseline active, the "ALL off +
% channel disabled + AWGN off" row may show a residual BER of order 1e-4
% rather than exactly 0. This is NOT a bug -- it is the unavoidable
% interaction between the CP-correlation flat-top (sometimes picks
% s_hat = +/- 1 instead of 0) and the baseline LS-with-complex-linear
% interpolation channel estimator at the band-edge SCs. Replacing either
% block (phase-aware interpolation, or a sharper STO estimator) drives it
% to exact zero. It does not affect any operating-condition BER.

clearvars; close all;

SNR_dB    = 28;             % single high-SNR point -- AWGN floor visible here
N_TRIALS  = 40;
seed0     = 1000;

% Each scenario: name, fd override, cfo override, sto override.
% Empty [] = use the default (random / 3000 Hz).
% Each scenario row: name, fd, cfo, sto, disable_channel, snr_override.
% Empty [] = use default (random / 3000 Hz / channel on / use SNR_dB arg).
scen = { ...
    'baseline (all impairments on)',                 [],   [],   [], false, []  ; ...
    'Doppler off only',                               0,   [],   [], false, []  ; ...
    'STO off only',                                  [],   [],    0, false, []  ; ...
    'CFO off only',                                  [],    0,   [], false, []  ; ...
    'STO + CFO off (Doppler still 3000)',            [],    0,    0, false, []  ; ...
    'ALL impairments off (fd=0, cfo=0, sto=0)',       0,    0,    0, false, []  ; ...
    'ALL off + channel disabled',                     0,    0,    0, true,  []  ; ...
    'ALL off + channel disabled + AWGN off',          0,    0,    0, true,  Inf };

fprintf('Diagnostic at Eb/N0 = %d dB (unless overridden), %d trials/scenario\n\n', ...
        SNR_dB, N_TRIALS);
fprintf('%-46s   %s\n', 'Scenario', 'BER');
fprintf('%-46s   %s\n', repmat('-', 1, 46), '----------');

for s = 1:size(scen, 1)
    name = scen{s, 1};
    opts = struct();
    if ~isempty(scen{s, 2}), opts.fd_doppler      = scen{s, 2}; end
    if ~isempty(scen{s, 3}), opts.cfo_Hz          = scen{s, 3}; end
    if ~isempty(scen{s, 4}), opts.sto_samples     = scen{s, 4}; end
    if scen{s, 5},           opts.disable_channel = true;       end
    if isempty(scen{s, 6}), snr_use = SNR_dB; else, snr_use = scen{s, 6}; end

    n_err = 0; n_bits = 0;
    for t = 1:N_TRIALS
        seed = seed0 + s*1e4 + t;
        [rx, params, ref] = ofdm_tx_and_channel(snr_use, seed, opts);
        bits = ofdm_rx(rx, params);
        L = min(length(bits), length(ref));
        n_err  = n_err  + sum(bits(1:L) ~= ref(1:L));
        n_bits = n_bits + L;
    end
    BER = n_err / max(n_bits, 1);
    fprintf('%-46s   %9.3e   (%d / %d)\n', name, BER, n_err, n_bits);
end
