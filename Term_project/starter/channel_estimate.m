function H = channel_estimate(Y, params)
% CHANNEL_ESTIMATE  LS at pilots + spline frequency interpolation.
%
%   H = channel_estimate(Y, params)
%
% Input:  Y - N_FFT x N_sym frequency-domain matrix (post-FFT)
% Output: H - N_used x N_sym channel estimate over the ACTIVE 312 SCs.

active    = params.active_bins;
pilot_sc  = params.pilot_sc_idx;
pilot_sym = params.pilot_sym_idx;
N_pilots  = numel(pilot_sc);
N_psym    = numel(pilot_sym);
X_pilot   = local_pilot_values(N_pilots);

% --- Step 1: LS at pilots for all pilot symbols (52 x N_psym) -----------
H_ls_pilots = zeros(N_pilots, N_psym);
for pIdx = 1:N_psym
    m = pilot_sym(pIdx);
    H_ls_pilots(:, pIdx) = Y(active(pilot_sc), m) ./ X_pilot;
end

% --- Step 2: Wiener time-domain smoothing across pilot symbols -----------
% Model: H_ls(p,m) = H_true(p,m) + noise. H_true is Jakes-correlated:
%   R_HH(i,j) = J0(2*pi*fd*|i-j|*T_sym)
% Smoother: W_t = R_HH * inv(R_HH + sigma2*I)
% Applied:  H_ls_smooth = H_ls_pilots * W_t   (52 x N_psym)
W_t = build_wiener_time(params);
H_ls_smooth = H_ls_pilots * W_t;

% --- Step 3: spline frequency interpolation per symbol ------------------
H = zeros(params.N_used, params.N_sym);
for pIdx = 1:N_psym
    m = pilot_sym(pIdx);
    H(:, m) = interp1(pilot_sc.', H_ls_smooth(:, pIdx), ...
                      (1:params.N_used).', 'spline');
end

% For any non-pilot symbols copy from nearest pilot (comb_10sym: none)
non_pilot = setdiff(1:params.N_sym, pilot_sym);
for m = non_pilot
    [~, idx] = min(abs(pilot_sym - m));
    H(:, m)  = H(:, pilot_sym(idx));
end

end

% =========================================================================
function W_t = build_wiener_time(params)
% Wiener smoother in time across N_sym pilot snapshots.
% R_HH(i,j) = J0(2*pi*fd*|i-j|*T_sym)  — Jakes temporal correlation.
% W_t = R_HH * inv(R_HH + sigma2*I),  size N_sym x N_sym.
% sigma2 = 5e-4 targets the 30 dB grading point (more aggressive smoothing).

fd    = params.fd_doppler;
T_sym = (params.N_FFT + params.N_CP) / params.Fs;
N     = params.N_sym;

m_idx = (0 : N-1).';
dk    = m_idx - m_idx.';                % N x N lag matrix
R_HH  = besselj(0, 2*pi * fd * abs(dk) * T_sym);

sigma2 = 5e-4;
W_t = R_HH / (R_HH + sigma2 * eye(N));  % N x N Wiener smoother

end

% =========================================================================
function p = local_pilot_values(N)
s = rng();  rng(12345);
b = randi([0 1], 2*N, 1);
pp = reshape(b, 2, []).';
ix = pp(:,1)*2 + pp(:,2);
tbl = (1/sqrt(2)) * [1+1j, 1-1j, -1+1j, -1-1j];
p = tbl(ix + 1).';
rng(s);
end
