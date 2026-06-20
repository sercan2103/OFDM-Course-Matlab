function [rx_samples, params, ref_bits] = ofdm_tx_and_channel(SNR_dB, seed, opts)
% OFDM_TX_AND_CHANNEL  Generate one received OFDM frame (black box).
%
%   [rx_samples, params, ref_bits] = ofdm_tx_and_channel(SNR_dB, seed)
%   [rx_samples, params, ref_bits] = ofdm_tx_and_channel(SNR_dB, seed, opts)
%
% Pipeline: bit-gen -> QPSK map -> pilot insertion -> IFFT -> CP add
%        -> FPV_C multipath (Jakes Doppler) -> STO injection
%        -> CFO injection -> AWGN.
%
% Inputs:
%   SNR_dB - Eb/N0 in dB
%   seed   - integer; seeds bits, channel, STO, and CFO together
%   opts   - (optional) struct with override fields for debugging:
%              .fd_doppler      - override maximum Doppler shift (Hz)
%              .cfo_Hz          - override CFO value (Hz, fixed instead of random)
%              .sto_samples     - override STO value (samples, fixed instead of random)
%              .disable_channel - true => bypass the FPV_C multipath block
%                                 entirely (signal passes through clean)
%            Any field that is missing or empty falls back to the defaults
%            from params_data() / the original randomised injection.
%
%   Passing SNR_dB = Inf disables AWGN as well (debug only). For the
%   official BER sweep the defaults are: FPV_C multipath ON, random STO in
%   [-32, +32] samples, random CFO in [-12, +12] kHz, AWGN at the swept
%   Eb/N0. Only fd_doppler is intended to be set during normal evaluation
%   (sweep over {0, 1000, 2000, 3000} Hz).
%
% Outputs:
%   rx_samples - complex baseband, column vector
%   params     - frame parameter struct (= params_data(), with overrides
%                applied to params.fd_doppler if requested)
%   ref_bits   - ground-truth bits, column vector (for BER scoring)

if nargin < 3, opts = struct(); end

rng(seed);
params = params_data();

% --- Apply optional Doppler override -------------------------------------
if isfield(opts, 'fd_doppler') && ~isempty(opts.fd_doppler)
    params.fd_doppler = opts.fd_doppler;
end

% --- 1) bit generation --------------------------------------------------
N_pilot_syms = numel(params.pilot_sym_idx);
N_data_per_pilot_sym = numel(params.data_sc_idx);
N_data_per_non_pilot = params.N_used;
N_qpsk = N_pilot_syms * N_data_per_pilot_sym + ...
         (params.N_sym - N_pilot_syms) * N_data_per_non_pilot;
N_bits = N_qpsk * 2;
ref_bits = randi([0 1], N_bits, 1);

% --- 2) QPSK mapping (gray, unit average power) -------------------------
pairs = reshape(ref_bits, 2, []).';
idx   = pairs(:,1)*2 + pairs(:,2);
qpsk_tbl = (1/sqrt(2)) * [1+1j, 1-1j, -1+1j, -1-1j];
qpsk_syms = qpsk_tbl(idx + 1).';

% --- 3) Pilot QPSK values (known, deterministic) ------------------------
pilot_vals = local_pilot_values(numel(params.pilot_sc_idx));

% --- 4) Build resource grid (N_FFT x N_sym) -----------------------------
grid = zeros(params.N_FFT, params.N_sym);
cur = 1;
for m = 1:params.N_sym
    if any(m == params.pilot_sym_idx)
        n_data = N_data_per_pilot_sym;
        d = qpsk_syms(cur:cur+n_data-1);  cur = cur + n_data;
        grid(params.active_bins(params.pilot_sc_idx), m) = pilot_vals;
        grid(params.active_bins(params.data_sc_idx),  m) = d;
    else
        n_data = N_data_per_non_pilot;
        d = qpsk_syms(cur:cur+n_data-1);  cur = cur + n_data;
        grid(params.active_bins, m) = d;
    end
end

% --- 5) IFFT + CP per symbol --------------------------------------------
tx_time = zeros((params.N_FFT + params.N_CP) * params.N_sym, 1);
for m = 1:params.N_sym
    s   = ifft(grid(:, m), params.N_FFT) * sqrt(params.N_FFT);
    cp  = s(end - params.N_CP + 1 : end);
    tx_time((m-1)*(params.N_FFT+params.N_CP)+1 : m*(params.N_FFT+params.N_CP)) = [cp; s];
end

% --- 6) Multipath channel (Rayleigh + Jakes Doppler) --------------------
% Debug switch: skip the channel entirely (signal passes through clean).
% comm.RayleighChannel requires MaximumDopplerShift > 0 to instantiate the
% Jakes spectrum. For fd = 0 we fall back to a static frequency-selective
% Rayleigh realisation (single draw of the tap gains, no time variation).
if isfield(opts, 'disable_channel') && opts.disable_channel
    ch_out = tx_time;
elseif params.fd_doppler > 0
    ch = comm.RayleighChannel( ...
        'SampleRate',            params.Fs, ...
        'PathDelays',            params.ch_delays_ns * 1e-9, ...
        'AveragePathGains',      params.ch_gains_dB, ...
        'MaximumDopplerShift',   params.fd_doppler, ...
        'DopplerSpectrum',       doppler('Jakes'), ...
        'RandomStream',          'mt19937ar with seed', ...
        'Seed',                  seed);
    ch_out = ch(tx_time);
else
    % Static multipath: convolve with one Rayleigh draw of each tap.
    delay_samps = round(params.ch_delays_ns * 1e-9 * params.Fs);
    lin_gains   = 10.^(params.ch_gains_dB / 20);
    taps        = zeros(max(delay_samps) + 1, 1);
    s_state     = rng(); rng(seed + 7919);   % decouple from bit stream RNG
    for k = 1:numel(delay_samps)
        h_k = (randn + 1j*randn) / sqrt(2);  % unit-variance complex Gaussian
        taps(delay_samps(k) + 1) = taps(delay_samps(k) + 1) + lin_gains(k) * h_k;
    end
    rng(s_state);
    ch_out = filter(taps, 1, tx_time);
end

% --- 7) STO injection (sample shift in [-32, +32]) ----------------------
if isfield(opts, 'sto_samples') && ~isempty(opts.sto_samples)
    sto = opts.sto_samples;
else
    sto = randi(params.sto_range_samples);
end
if sto >= 0
    rx_shifted = [zeros(sto, 1); ch_out];
else
    rx_shifted = ch_out(-sto + 1 : end);
end

% --- 8) CFO injection (uniform in [-12 kHz, +12 kHz]) -------------------
if isfield(opts, 'cfo_Hz') && ~isempty(opts.cfo_Hz)
    cfo_Hz = opts.cfo_Hz;
else
    cfo_Hz = (rand*2 - 1) * params.cfo_range_Hz(2);
end
n_idx  = (0:length(rx_shifted)-1).';
rx_cfo = rx_shifted .* exp(1j * 2 * pi * cfo_Hz * n_idx / params.Fs);

% --- 9) AWGN at the requested Eb/N0 -------------------------------------
% QPSK: Es = 2*Eb.  Scale for the fact that only N_used / N_FFT of the
% available bandwidth carries energy. SNR_dB = Inf disables AWGN (debug).
if isinf(SNR_dB)
    rx_samples = rx_cfo;
else
    EsN0_dB    = SNR_dB + 10*log10(2);
    SNR_eff_dB = EsN0_dB + 10*log10(params.N_used / params.N_FFT);
    rx_samples = awgn(rx_cfo, SNR_eff_dB, 'measured');
end

end

% =========================================================================
function p = local_pilot_values(N)
% Deterministic QPSK pilot sequence (same every frame; known to the RX).
s = rng();  rng(12345);
b = randi([0 1], 2*N, 1);
pp = reshape(b, 2, []).';
ix = pp(:,1)*2 + pp(:,2);
tbl = (1/sqrt(2)) * [1+1j, 1-1j, -1+1j, -1-1j];
p = tbl(ix + 1).';
rng(s);
end
