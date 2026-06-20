clear; clc; close all;

N      = 64;
Ncp    = 16;
h      = [1 0 0 0.6];
h_long = [ones(1,5) zeros(1,5) 0.8*ones(1,5) zeros(1,5) 0.5];
Nbits  = 2e5;
rng(42);

H     = fft(h, N);
Hlong = fft(h_long, N);

% =========================================================================
% TODO 1
% =========================================================================
bits = randi([0 1], Nbits, 1);

snr1     = 0:2:12;
snr2     = 0:2:30;
ber_awgn = zeros(size(snr1));
ber_nocp = zeros(size(snr2));

for i = 1:length(snr1)
    ber_awgn(i) = ofdm_chain(bits, snr1(i), 1, N, Ncp, ones(N,1), true);
end
for i = 1:length(snr2)
    ber_nocp(i) = ofdm_chain(bits, snr2(i), h, N, 0, H.', true);
end

snr_th = linspace(0,12,300);
ber_th = 0.5*erfc(sqrt(10.^(snr_th/10)));

bits_c     = randi([0 1], 8*N*2, 1);
con_nocp   = ofdm_chain_constel(bits_c, 20, h, N, 0,   H.', false);

figure('Color','w','Position',[50 50 1100 460]);

subplot(1,2,1);
semilogy(snr_th, ber_th, 'b-', 'LineWidth',2, 'DisplayName','QPSK theory (AWGN)'); hold on;
semilogy(snr1, ber_awgn, 'ko', 'MarkerFaceColor','k', 'MarkerSize',7, ...
    'DisplayName','Phase 2: AWGN sanity (CP on)');
semilogy(snr2, ber_nocp, 'rs-', 'MarkerFaceColor','r', 'MarkerSize',7, 'LineWidth',1.8, ...
    'DisplayName','Phase 3: multipath, NO CP (the shock)');
grid on; grid minor;
xlim([0 30]); ylim([1e-6 1]);
xlabel('E_b/N_0 [dB]'); ylabel('BER');
title('Build (AWGN) \rightarrow Break (multipath, no CP)');
legend('Location','southwest','FontSize',9);

subplot(1,2,2);
s = 1/sqrt(2);
plot(real(con_nocp), imag(con_nocp), '.', 'Color',[0.9 0.3 0.3], 'MarkerSize',3); hold on;
plot(real(s*[1+1j;1-1j;-1+1j;-1-1j]), imag(s*[1+1j;1-1j;-1+1j;-1-1j]), ...
    'rx', 'MarkerSize',16, 'LineWidth',2.5);
axis equal; grid on;
xlim([-2 2]); ylim([-2 2]);
xlabel('Re'); ylabel('Im');
title('Constellation @ 20 dB (no CP)');
text(0,-1.8,'BER \approx 10^{-2} — Chain UNUSABLE', ...
    'HorizontalAlignment','center','Color','r','FontWeight','bold');

sgtitle('TODO 1 — Build the Chain & Witness the BER Floor','FontWeight','bold');

% =========================================================================
% TODO 2
% =========================================================================
snr_t2  = 0:1:20;
ber_off = zeros(size(snr_t2));
ber_on  = zeros(size(snr_t2));

for i = 1:length(snr_t2)
    ber_off(i) = ofdm_chain(bits, snr_t2(i), h, N, 0,   H.', true);
    ber_on(i)  = ofdm_chain(bits, snr_t2(i), h, N, Ncp, H.', true);
end

snr_th2    = linspace(0,20,300);
ber_th2    = 0.5*erfc(sqrt(10.^(snr_th2/10)));
tax_dB     = 10*log10((N+Ncp)/N);
ber_th2_cp = 0.5*erfc(sqrt(10.^((snr_th2-tax_dB)/10)));

bits_c2 = randi([0 1], 40*N*2, 1);
con_on  = ofdm_chain_constel(bits_c2, 15, h, N, Ncp, H.', true);

figure('Color','w','Position',[50 50 1100 500]);

subplot(1,2,1);
semilogy(snr_th2, ber_th2,    'k-',  'LineWidth',1.5, 'DisplayName','AWGN theory (reference)'); hold on;
semilogy(snr_th2, ber_th2_cp, 'k--', 'LineWidth',1.2, ...
    'DisplayName',sprintf('Theory + %.1f dB CP tax', tax_dB));
semilogy(snr_t2, max(ber_off,1e-6), 'rs-', 'MarkerFaceColor','r', 'MarkerSize',7, ...
    'LineWidth',1.8, 'DisplayName','from TODO 1: CP OFF');
semilogy(snr_t2, max(ber_on,1e-6),  'go-', 'MarkerFaceColor','g', 'MarkerSize',7, ...
    'LineWidth',1.8, 'DisplayName','TODO 2: CP ON (N_{cp}=16)');
grid on; grid minor;
xlim([0 20]); ylim([1e-6 1]);
xlabel('E_b/N_0 [dB]'); ylabel('BER');
title('TODO 2 — One Switch, BER Snaps onto Theory');
legend('Location','southwest','FontSize',9);
text(6, 3e-4, '\leftarrow ~1 dB CP tax','Color',[0.4 0.4 0.4]);

subplot(1,2,2);
plot(real(con_on), imag(con_on), '.', 'Color',[0.1 0.7 0.2], 'MarkerSize',3); hold on;
viscircles([0 0], 1.0, 'Color',[0.7 0.7 0.7], 'LineWidth',0.5);
axis equal; grid on;
xlim([-1.8 1.8]); ylim([-1.8 1.8]);
xlabel('Re'); ylabel('Im');
title('Constellation @ 15 dB (CP ON, EQ ON)');
text(0,-1.6,'Clean QPSK — BER < 10^{-4}', ...
    'HorizontalAlignment','center','Color',[0 0.5 0],'FontWeight','bold');

sgtitle('TODO 2 — Flip CP On: BER Restores to Theory','FontWeight','bold');
fprintf('CP tax = %.2f dB\n', tax_dB);

% =========================================================================
% TODO 3
% =========================================================================
ncp_vec   = [0 4 8 12 16 20 24 28 32];
tau_long  = length(h_long) - 1;
ber_sweep = zeros(size(ncp_vec));
bits_long = randi([0 1], 5e5, 1);

for i = 1:length(ncp_vec)
    ber_sweep(i) = ofdm_chain(bits_long, 10, h_long, N, ncp_vec(i), Hlong.', true);
end

figure('Color','w','Position',[100 100 750 500]);
semilogy(ncp_vec, max(ber_sweep,1e-6), 'o-', 'Color',[0.85 0.45 0], ...
    'MarkerFaceColor',[0.85 0.45 0], 'MarkerSize',9, 'LineWidth',2.2, ...
    'DisplayName','BER'); hold on;
xline(tau_long,'b--','LineWidth',2,'DisplayName',sprintf('N_{cp} = \\tau = %d',tau_long));
fill([0 tau_long tau_long 0],[1e-6 1e-6 1 1],[1 0.85 0.85],'FaceAlpha',0.2,'EdgeColor','none');
fill([tau_long 33 33 tau_long],[1e-6 1e-6 1 1],[0.85 1 0.85],'FaceAlpha',0.2,'EdgeColor','none');
text(tau_long/2, 0.4,'ISI / ICI','HorizontalAlignment','center','Color',[0.8 0.1 0.1],'FontWeight','bold');
text(tau_long+5, 2e-5,'safe','HorizontalAlignment','center','Color',[0 0.55 0],'FontWeight','bold');
text(tau_long+0.4, 0.25,sprintf('knee at\nN_{cp}=\\tau=%d',tau_long),'Color','b','FontWeight','bold');
grid on; grid minor;
xlim([-1 33]); ylim([1e-6 1]);
xlabel('CP length N_{cp} [samples]'); ylabel('BER  (E_b/N_0 = 10 dB)');
title(sprintf('TODO 3 — BER vs N_{cp}  (h_{long},  \\tau = %d)', tau_long));
legend('Location','northeast');

% =========================================================================
% What-if: remove EQ
% =========================================================================
bits_wi = randi([0 1], 60*N*2, 1);
Ymat    = ofdm_chain_Ymat(bits_wi, 15, h, N, Ncp);

figure('Color','w','Position',[100 100 650 580]);
cols = lines(8);
hold on;
for k = 1:8
    plot(real(Ymat(k,:)), imag(Ymat(k,:)), '.', 'Color',cols(k,:), ...
        'MarkerSize',5, 'DisplayName',sprintf('k = %d',k));
end
plot(real(s*[1+1j;1-1j;-1+1j;-1-1j]), imag(s*[1+1j;1-1j;-1+1j;-1-1j]), ...
    'kx', 'MarkerSize',14,'LineWidth',2.5,'HandleVisibility','off');
viscircles([0 0],1.0,'Color',[0.7 0.7 0.7],'LineWidth',0.4);
axis equal; grid on;
xlim([-1.8 1.8]); ylim([-1.8 1.8]);
xlabel('Re'); ylabel('Im');
title({'What-If — Y[k] with NO Equalizer', ...
       '(each subcarrier rotated + scaled by its own H[k])'});
legend('Location','bestoutside','NumColumns',2,'FontSize',9);

% =========================================================================
% Local functions
% =========================================================================
function syms = qpsk_map(bits)
    bits = bits(:);
    p    = reshape(bits, 2, []).';
    syms = ((1 - 2*double(p(:,1))) + 1j*(1 - 2*double(p(:,2)))) / sqrt(2);
end

function bits = qpsk_demap(syms)
    syms = syms(:);
    bits = reshape([double(real(syms)<0) double(imag(syms)<0)].', [], 1);
end

function ber = ofdm_chain(bits, EbN0_dB, hch, N, ncp, Hk, eq_on)
    syms = qpsk_map(bits);
    Nblk = floor(length(syms)/N);
    syms = syms(1:Nblk*N);
    bits = bits(1:Nblk*N*2);

    tx = ifft(reshape(syms, N, Nblk), N) * sqrt(N);
    if ncp > 0
        tx = [tx(end-ncp+1:end,:); tx];
    end

    rx = conv(tx(:), hch);
    rx = rx(1:numel(tx));

    snr = 10^(EbN0_dB/10) * 2 * N / (N + ncp);
    rx  = rx + sqrt(1/(2*snr)) * (randn(size(rx)) + 1j*randn(size(rx)));

    rx = reshape(rx, N+ncp, Nblk);
    if ncp > 0
        rx = rx(ncp+1:end,:);
    end

    Y = fft(rx, N) / sqrt(N);
    if eq_on
        Y = Y ./ repmat(Hk(:), 1, Nblk);
    end

    rx_bits = qpsk_demap(Y(:));
    Nb  = min(length(bits), length(rx_bits));
    ber = sum(bits(1:Nb) ~= rx_bits(1:Nb)) / Nb;
end

function constel = ofdm_chain_constel(bits, EbN0_dB, hch, N, ncp, Hk, eq_on)
    syms = qpsk_map(bits);
    Nblk = floor(length(syms)/N);

    tx = ifft(reshape(syms(1:Nblk*N), N, Nblk), N) * sqrt(N);
    if ncp > 0
        tx = [tx(end-ncp+1:end,:); tx];
    end

    rx = conv(tx(:), hch);
    rx = rx(1:numel(tx));

    snr = 10^(EbN0_dB/10) * 2 * N / (N + ncp);
    rx  = rx + sqrt(1/(2*snr)) * (randn(size(rx)) + 1j*randn(size(rx)));

    rx = reshape(rx, N+ncp, Nblk);
    if ncp > 0
        rx = rx(ncp+1:end,:);
    end

    Y = fft(rx, N) / sqrt(N);
    if eq_on
        Y = Y ./ repmat(Hk(:), 1, Nblk);
    end
    constel = Y(:);
end

function Ymat = ofdm_chain_Ymat(bits, EbN0_dB, hch, N, ncp)
    syms = qpsk_map(bits);
    Nblk = floor(length(syms)/N);

    tx = ifft(reshape(syms(1:Nblk*N), N, Nblk), N) * sqrt(N);
    if ncp > 0
        tx = [tx(end-ncp+1:end,:); tx];
    end

    rx = conv(tx(:), hch);
    rx = rx(1:numel(tx));

    snr = 10^(EbN0_dB/10) * 2 * N / (N + ncp);
    rx  = rx + sqrt(1/(2*snr)) * (randn(size(rx)) + 1j*randn(size(rx)));

    rx = reshape(rx, N+ncp, Nblk);
    if ncp > 0
        rx = rx(ncp+1:end,:);
    end

    Ymat = fft(rx, N) / sqrt(N);
end
