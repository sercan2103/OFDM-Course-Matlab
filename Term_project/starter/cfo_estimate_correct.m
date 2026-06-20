function rx_out = cfo_estimate_correct(rx_in, params)
% CFO_ESTIMATE_CORRECT  Residual CFO estimation via CP self-correlation.
%
%   rx_out = cfo_estimate_correct(rx_in, params)
%
% Method: van de Beek 1997 CP self-correlation.
%   The CP samples r[n] match r[n + N_FFT] up to a phase rotation of
%   exp(-j*2*pi*cfo*N_FFT/Fs). Averaging across all 10 symbols and
%   taking the angle gives an unbiased estimate:
%
%     cfo_hat = -angle( sum_{m,n} r[n] * conj(r[n+N_FFT]) ) * Fs / (2*pi*N_FFT)
%
%   Unambiguous range: +/- Fs/(2*N_FFT) = +/- 15 kHz > spec +/- 12 kHz.
%   Unbiased under multipath: all FPV_C taps (max delay 12 samples) are
%   well within the 288-sample CP, so the CP-body match is exact per tap.
%
% NOTE: A two-stage pilot-cross-correlation refinement was tested and
% reverted. At high SNR the van de Beek residual is already < 50 Hz;
% the pilot refinement adds Doppler-induced channel phase noise (~200 Hz
% std at fd=2000 Hz) that exceeds the residual it is trying to correct.

N_FFT   = params.N_FFT;
N_CP    = params.N_CP;
Fs      = params.Fs;
sym_len = N_FFT + N_CP;

% Average CP self-correlation across all available symbols
n_syms = floor(length(rx_in) / sym_len);
c      = complex(0);
for m = 1 : n_syms
    base = (m - 1) * sym_len;
    if base + N_FFT + N_CP > length(rx_in), break; end
    c = c + sum(rx_in(base + 1 : base + N_CP) .* ...
                conj(rx_in(base + N_FFT + 1 : base + N_FFT + N_CP)));
end

cfo_hat = -angle(c) * Fs / (2 * pi * N_FFT);

% Derotate sample-by-sample
n      = (0 : length(rx_in) - 1).';
rx_out = rx_in .* exp(-1j * 2 * pi * cfo_hat * n / Fs);

end
