function bits = qpsk_demap(X_hat, params)
% QPSK_DEMAP  Hard-decision QPSK demapping (Gray-coded).
%
%   bits = qpsk_demap(X_hat, params)
%
% Input:  X_hat - N_used x N_sym equalised symbols on the 312 active SCs.
% Output: bits  - column vector, in the same order ofdm_tx_and_channel
%                 generated them.
%
% On pilot-bearing OFDM symbols (1, 4, 7, 10) only the 156 data SCs carry
% bits; the 156 pilot SCs are skipped here. On non-pilot OFDM symbols all
% 312 SCs are data.

bits = [];
for m = 1:params.N_sym
    if any(m == params.pilot_sym_idx)
        d = X_hat(params.data_sc_idx, m);
    else
        d = X_hat(:, m);
    end
    % Gray-coded QPSK: bit_msb = sign(real)<0; bit_lsb = sign(imag)<0
    b_msb = real(d) < 0;
    b_lsb = imag(d) < 0;
    bb = [b_msb, b_lsb].';
    bits = [bits; bb(:)];
end

bits = double(bits);

end
