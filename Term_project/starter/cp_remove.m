function rx_no_cp = cp_remove(rx_in, params)
% CP_REMOVE  Strip the CP from each OFDM symbol.
%
%   rx_no_cp = cp_remove(rx_in, params)
%
% Output: N_FFT x N_sym matrix of complex time-domain samples (no CP).
%
% (You may modify this, but the standard interpretation -- skip the first
%  N_CP samples of each (N_CP + N_FFT) block -- is essentially the only
%  thing that makes sense if no other corrections were applied.)

N_FFT   = params.N_FFT;
N_CP    = params.N_CP;
N_sym   = params.N_sym;
sym_len = N_FFT + N_CP;

% Truncate or pad rx_in so it holds exactly N_sym symbol blocks
needed = sym_len * N_sym;
if length(rx_in) < needed
    rx_in(end+1:needed) = 0;
else
    rx_in = rx_in(1:needed);
end

rx_no_cp = zeros(N_FFT, N_sym);
for m = 1:N_sym
    blk = rx_in((m-1)*sym_len + 1 : m*sym_len);
    rx_no_cp(:, m) = blk(N_CP + 1 : end);    % skip CP
end

end
