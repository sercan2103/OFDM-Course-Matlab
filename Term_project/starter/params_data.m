function P = params_data()
% PARAMS_DATA  Frame parameters for the OFDM term project.
%
% Returns a struct used by both the transmitter (ofdm_tx_and_channel) and
% the receiver (ofdm_rx) and all sub-functions. You do not need to change
% this file. It is the single source of truth for frame geometry.

% --- Core OFDM numerology -------------------------------------------------
P.Fs        = 30.72e6;           % sampling frequency (Hz)
P.N_FFT     = 1024;              % FFT size (samples)
P.N_CP      = 288;               % CP length (samples)
P.N_sym     = 10;                % OFDM symbols per frame
P.SCS       = 30e3;              % sub-carrier spacing (Hz)
P.BW        = 10e6;              % occupied bandwidth (Hz)

% --- Resource blocks ------------------------------------------------------
P.N_RB      = 26;                % number of resource blocks
P.SC_per_RB = 12;                % sub-carriers per resource block
P.N_used    = P.N_RB * P.SC_per_RB;   % 312 active SCs (DC excluded)

% --- Active sub-carrier mapping ------------------------------------------
% FFT-bin indices of the 312 active SCs in MATLAB's natural FFT order.
% DC = bin 1.  Positive freqs: bins 2..157.  Negative freqs: bins 869..1024.
half = P.N_used / 2;                                              % 156
P.active_bins = [2:(half+1), (P.N_FFT - half + 1):P.N_FFT];       % length 312

% --- Pilot pattern: comb_10sym_64 ----------------------------------------
% Interleaved 3+2 pilot pattern hitting the 64-pilot cap exactly:
%   Even RBs 0,2,...,22 (12 RBs): pilots at {1,5,9}  — spacing 4 SCs
%   Odd  RBs 1,3,...,23 + RBs 24,25 (14 RBs): pilots at {1,7} — spacing 6 SCs
%   Total: 12×3 + 14×2 = 64 pilot SCs / symbol (= cap)
%
% Average pilot spacing = 312/64 ≈ 4.9 SCs vs 6 SCs before → better
% frequency interpolation. All 10 symbols carry pilots (comb pattern kept).
%
% TX modification note: extended from comb_10sym (52 pilots, spacing 6) to
% hit the 64-pilot cap using an interleaved 3+2 RB layout.
P.pilot_offsets_in_RB = [1, 7];       % default 2-pilot offsets (reference only)
P.pilot_sym_idx       = 1 : P.N_sym;  % all 10 symbols

pilot_pos = [];
for rb = 0:(P.N_RB - 1)
    if mod(rb, 2) == 0 && rb <= 22    % even RBs 0..22: 3 pilots
        offsets = [1, 5, 9];
    else                               % odd RBs + RBs 24,25: 2 pilots
        offsets = [1, 7];
    end
    pilot_pos = [pilot_pos, rb * P.SC_per_RB + offsets];
end
P.pilot_sc_idx = pilot_pos;           % length = 64
P.data_sc_idx  = setdiff(1:P.N_used, P.pilot_sc_idx);  % length = 248

% --- Channel (FPV_C urban / suburban) ------------------------------------
P.ch_delays_ns = [0,  50, 150,  400];
P.ch_gains_dB  = [0,  -5, -10,  -16];
P.fd_doppler   = 3000;            % maximum Doppler shift (Hz, Jakes)
P.fc           = 2.4e9;           % carrier frequency

% --- Impairment ranges (instructor picks the actual values) --------------
P.sto_range_samples = [-32, 32];  % residual STO bound (samples)
P.cfo_range_Hz      = [-12e3, 12e3];  % residual CFO bound (Hz)

end
