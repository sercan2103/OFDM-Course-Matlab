function Y = ofdm_demod(rx_no_cp, params)
% OFDM_DEMOD  Per-symbol FFT.
%
%   Y = ofdm_demod(rx_no_cp, params)
%
% Input:  N_FFT x N_sym time-domain matrix (CP already removed)
% Output: N_FFT x N_sym frequency-domain matrix (DC at bin 1)
%
% Energy is normalised so unit-power TX symbols yield unit-power FFT bins.

Y = fft(rx_no_cp, params.N_FFT, 1) / sqrt(params.N_FFT);

end
