function X_hat = equalize(Y, H, params)
% EQUALIZE  ZF + decision-directed ICI cancellation.
%
%   X_hat = equalize(Y, H, params)
%
% ICI model: for a linearly time-varying channel within one OFDM symbol,
% the received SC-k signal is
%
%   Y[k] = H_avg[k]*X[k]
%         + sum_{Dk!=0} dH[k-Dk] * w(Dk) * X[k-Dk]  <- ICI
%         + N[k]
%
% where dH[l] = channel slope at SC l within the symbol (estimated from
% adjacent symbol channel estimates), and
%
%   w(Dk) = 1 / (N_FFT * (exp(-j*2*pi*Dk/N_FFT) - 1))
%
% derived by DFT-ing the linear ramp (n/N - 1/2) over n=0..N-1.
%
% Algorithm:
%   1) ZF equalize  ->  X_zf
%   2) Hard QPSK decisions  ->  X_dec
%   3) Estimate dH per symbol from adjacent H snapshots
%   4) Subtract ICI from Y (truncated to ±Q SCs)
%   5) Re-equalize  ->  X_hat

Y_active = Y(params.active_bins, :);   % N_used x N_sym
N_used   = params.N_used;
N_sym    = params.N_sym;
N_FFT    = params.N_FFT;
Q        = 6;   % ICI cancellation radius (SCs each side)

% --- Precompute ICI weights w(Dk) for Dk = +1..+Q and -1..-Q -----------
dk    = (1:Q);
w_pos = 1 ./ (N_FFT * (exp(-1j*2*pi*dk/N_FFT) - 1));  % w(+1..+Q)
w_neg = 1 ./ (N_FFT * (exp( 1j*2*pi*dk/N_FFT) - 1));  % w(-1..-Q)

% --- Stage 1: initial ZF -------------------------------------------------
X_zf = Y_active ./ H;

% --- Stage 2: hard QPSK decisions ----------------------------------------
s     = 1/sqrt(2);
X_dec = s * (sign(real(X_zf)) + 1j*sign(imag(X_zf)));  % N_used x N_sym

% --- Stage 3: channel slope estimate per symbol --------------------------
% dH(:,m) approximates the total channel change across symbol m's useful
% duration, estimated from central difference of adjacent H snapshots.
% Scale by T_useful/T_symbol = N_FFT/(N_FFT+N_CP) to convert from
% symbol-to-symbol difference to within-symbol variation.
scale = params.N_FFT / (params.N_FFT + params.N_CP);  % ~0.78

dH = zeros(N_used, N_sym);
dH(:, 1)           = (H(:, 2)     - H(:, 1))           * scale;
dH(:, N_sym)       = (H(:, N_sym) - H(:, N_sym-1))     * scale;
dH(:, 2:N_sym-1)   = (H(:, 3:N_sym) - H(:, 1:N_sym-2)) * (scale/2);

% --- Stage 4: compute and subtract ICI -----------------------------------
Y_corr = Y_active;
for m = 1:N_sym
    ici   = zeros(N_used, 1);
    dHm   = dH(:, m);
    Xdm   = X_dec(:, m);

    for qi = 1:Q
        % ICI at dst from high-side src (src = dst+qi, Dk = -qi)
        src = (qi+1 : N_used);
        dst = src - qi;
        ici(dst) = ici(dst) + dHm(src) .* w_neg(qi) .* Xdm(src);

        % ICI at dst from low-side src (src = dst-qi, Dk = +qi)
        src = (1 : N_used-qi);
        dst = src + qi;
        ici(dst) = ici(dst) + dHm(src) .* w_pos(qi) .* Xdm(src);
    end
    Y_corr(:, m) = Y_active(:, m) - ici;
end

% --- Stage 5: final ZF on ICI-corrected Y --------------------------------
X_hat = Y_corr ./ H;

end
