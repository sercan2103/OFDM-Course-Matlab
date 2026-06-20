function bits = ofdm_rx(rx_samples, params)
% OFDM_RX  Top-level student receiver. Orchestrates the chain.
%
%   bits = ofdm_rx(rx_samples, params)
%
% The chain calls the following sub-functions, each of which lives in its
% own file. Some are provided as working starters (you may replace them
% for better performance); others are STUBS you must implement.
%
%   rx_samples  -> sto_estimate_correct      [STUB - implement]
%               -> cfo_estimate_correct      [STUB - implement]
%               -> cp_remove                 [given]
%               -> ofdm_demod (FFT)          [given]
%               -> channel_estimate          [given - replaceable]
%               -> doppler_track             [STUB - implement]
%               -> equalize                  [given - replaceable]
%               -> qpsk_demap                [given]
%               -> bits

% 1) Time-domain corrections
rx1 = sto_estimate_correct(rx_samples, params);
rx2 = cfo_estimate_correct(rx1, params);

% 2) Remove CP, FFT
rx_no_cp = cp_remove(rx2, params);
Y        = ofdm_demod(rx_no_cp, params);

% 3) Channel estimation + Doppler tracking
H_init = channel_estimate(Y, params);
H      = doppler_track(H_init, Y, params);

% 4) Equalize + demap
X_hat  = equalize(Y, H, params);
bits   = qpsk_demap(X_hat, params);

end
