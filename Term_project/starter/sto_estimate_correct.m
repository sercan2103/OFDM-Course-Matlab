function rx_out = sto_estimate_correct(rx_in, params)
% STO_ESTIMATE_CORRECT  Residual STO estimation via frequency-domain
%                       differential phase on pilot symbols.
%
%   rx_out = sto_estimate_correct(rx_in, params)
%
% A residual STO of delta samples causes a linear phase ramp across bins:
%
%   H_raw[k] = H_true[k] * exp(-j*2*pi*(k-1)*delta/N_FFT)
%
% For two pilot SCs with true frequency spacing dk (in SCS units):
%
%   dphi = angle(H_raw[k+dk] * conj(H_raw[k]))
%        ≈ -2*pi * dk * delta / N_FFT
%   =>  delta = -dphi * N_FFT / (2*pi * dk)
%
% Works for any pilot spacing (4 or 6 SCs in the mixed 3+2 RB pattern).
% Pairs crossing the DC guard band (|dk_true| ≈ 307) are detected via
% DC-centred frequency coordinates and skipped.
% Unambiguous range: +/- N_FFT / (2*dk_min) = +/- 128 samples >> spec +/- 32.
% CFO does not bias: constant phase per symbol cancels in within-symbol diff.

N_FFT    = params.N_FFT;
N_CP     = params.N_CP;
sym_len  = N_FFT + N_CP;
active   = params.active_bins;
pilot_sc = params.pilot_sc_idx;

% Known pilot QPSK values (must match the TX seed)
X_pilot = local_pilot_values(numel(pilot_sc));

% -------------------------------------------------------------------------
% DC-centred frequency positions of each pilot SC (in units of SCS).
%   Active SCs 1..156  → freq +1..+156  (positive half)
%   Active SCs 157..312 → freq -156..-1 (negative half)
% Consecutive pilot pairs that cross the DC gap appear to have index
% spacing ~6 but true frequency separation ~307 SCS — must be skipped.
% Threshold of 20 SCS cleanly separates in-band pairs (spacing 4 or 6)
% from the DC-crossing pair (spacing 307).
% -------------------------------------------------------------------------
half      = params.N_used / 2;                        % 156
sc_freq   = [(1:half), (-half:-1)];                   % 1x312
pf        = sc_freq(pilot_sc);                        % freq of each pilot SC
N_p       = numel(pilot_sc);

% Pre-select valid (non-DC-crossing) consecutive pairs
valid_pairs = false(1, N_p-1);
dk_pairs    = zeros(1, N_p-1);
for p = 1 : N_p-1
    dk = pf(p+1) - pf(p);
    if abs(dk) <= 20           % in-band pair (spacing 4 or 6)
        valid_pairs(p) = true;
        dk_pairs(p)    = dk;
    end
end

% Accumulate delta estimates across all pilot symbols
delta_accum = 0;
valid_syms  = 0;

for ps = 1 : numel(params.pilot_sym_idx)
    m       = params.pilot_sym_idx(ps);
    start_s = (m - 1) * sym_len + N_CP + 1;
    end_s   = start_s + N_FFT - 1;
    if end_s > length(rx_in), continue; end

    Y     = fft(rx_in(start_s : end_s), N_FFT) / sqrt(N_FFT);
    H_raw = Y(active(pilot_sc)) ./ X_pilot;

    % Weighted-average delta: each pair contributes its own spacing-normalised
    % differential phase estimate — handles mixed 4-SC and 6-SC spacings.
    delta_sum  = 0;
    pair_count = 0;
    for p = 1 : N_p-1
        if ~valid_pairs(p), continue; end
        dphi       = angle(H_raw(p+1) * conj(H_raw(p)));
        delta_sum  = delta_sum + (-dphi * N_FFT / (2*pi * dk_pairs(p)));
        pair_count = pair_count + 1;
    end

    delta_accum = delta_accum + delta_sum / max(pair_count, 1);
    valid_syms  = valid_syms + 1;
end

% Final integer estimate clamped to the spec'd range
delta_hat = round(delta_accum / max(valid_syms, 1));
delta_hat = max(min(delta_hat, 32), -32);

% -------------------------------------------------------------------------
% Apply integer sample shift to align the FFT window
% -------------------------------------------------------------------------
if delta_hat > 0
    rx_out = rx_in(delta_hat + 1 : end);           % skip delta_hat leading samples
elseif delta_hat < 0
    rx_out = [zeros(-delta_hat, 1); rx_in];         % prepend zeros
else
    rx_out = rx_in;
end

end

% =========================================================================
function p = local_pilot_values(N)
% Deterministic pilot sequence — must match ofdm_tx_and_channel.m exactly.
s = rng();  rng(12345);
b = randi([0 1], 2*N, 1);
pp = reshape(b, 2, []).';
ix = pp(:,1)*2 + pp(:,2);
tbl = (1/sqrt(2)) * [1+1j, 1-1j, -1+1j, -1-1j];
p = tbl(ix + 1).';
rng(s);
end
