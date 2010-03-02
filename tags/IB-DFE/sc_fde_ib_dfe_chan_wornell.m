% File: sc_fde_ib_dfe_chan_wornell.m
% ----------------------------------
% To compare with 'sc_fde_ib_dfe.m', this script implements SC-FDE-IB-DFE 
% system based on 'sc_fde_ib_dfe.m' and [1][2].
close all
clear all

%%
BW = 2e6; % bandwidth in bits/sec.
Ts = 1/BW; % symbol interval
% path_delays_norm = [0 2 5 7]; % multipath delay spread normalized to Ts
doppler_max = 0; % max doppler shift in Hz
% pdp = [0 0 -5 -5]; % power delay profile
path_delays_norm = [0: 1: 16];
pdp = zeros(1,length(path_delays_norm)); % severe-ISI, equally decaying Rayleigh fading channel

nfft = 128; % #of symbols in a fde block
length_cp = max(path_delays_norm); % #of symbols in cyclic prefix
length_pload = nfft-length_cp; % #of play load symbols in null UW scheme
nblock = 1; % #of fde blocks in each frame, >1 consider IBI; <=1 without considering IBI
numTrials = 1000; % #of frames in each evaluation
M_mod = 4; % Alphabet size
%%
% Compute error rate for different values of SNR.
SNR = 6:2:16; % Range of SNR values, in dB. receiver snr
numerr = zeros(4,length(SNR)); % zf, mmse, dfe, ib
ber_mfb = zeros(1,length(SNR));
h_waitbar = waitbar(0,'Please wait...');
for n = 1:length(SNR)
    chan = rayleighchan(Ts,doppler_max,path_delays_norm*Ts,pdp); % rayleigh fading channel
    chan.StorePathGains = true; % the complex path gain vector is stored 
                                % as the channel filter function processes 
                                % the signal
    chan.ResetBeforeFiltering = true; % each call to filter resets the 
                                      % state of chan before filtering
    chan.NormalizePathGains = true;  % If 1, the Rayleigh fading process 
                                      % is normalized such that the expected 
                                      % value of the path gains' total power is 1.
%     if (doppler_max > 0)
%         chan.StoreHistory = true; % to plot channel response
%     else
%         chan.StoreHistory = false;
%     end
%%
    for m = 1:numTrials % quasi-static channel modeling
        % Create a random digital message
        % Use M-QAM modulation.
%         x = randint(nfft,nblock,M_mod);
        x = randi([0,M_mod-1], length_pload, 1); % for UW scheme
        xtx_pload = modulate(modem.qammod(M_mod),x);
        xtxPower = sum(abs(xtx_pload).^2,1)/size(xtx_pload,1); % check transmitted signal power along each tx antenna
        xtx_pload = xtx_pload*diag(1./sqrt(xtxPower))./sqrt(1); % normalized and fixed total power
%         xtx = [xtx_pload(end-length_cp+1: end, :); xtx_pload]; % insert CP
        xtx = [xtx_pload; zeros(length_cp,1)]; % add UW at the end of a block.
                                          % the UW at the beginning of a
                                          % block is omitted, since the UW
                                          % we choosed here is a null
                                          % sequence, it will not interfere
                                          % with the pay load.
        XTX = fft(xtx,nfft); % DFT of transmitted symbols, used at rx end for filter design.
        M_XTX = XTX'*XTX/nfft; % Power of transmitted symbols, used at rx end for filter design.
%         xtx = reshape(xtx, [], 1);
        fadesig = filter(chan, xtx); % Effect of channel, quasi-static
        %%
        xtxPower = sum(abs(xtx(:)).^2)/length(xtx(:)); % check transmitted signal power
        xtxPower_dB = 10*log10(xtxPower); % convert to dB
        fadesigPower = sum(abs(fadesig(:)).^2)/length(fadesig(:)); % check fading signal power
        fadesigPower_dB = 10*log10(fadesigPower); % convert to dB
%         noisePower_dB = fadesigPower_dB-SNR(n); % compute noise power, in dB
        noisePower_dB = zeros(size(fadesigPower_dB))-SNR(n); % compute noise power, with 
                                                  % the assumption that the
                                                  % received signal at each
                                                  % rx-antenna is unit
                                                  % power (0dB).
        noisePower = 10^(noisePower_dB/10); % convert backto linear scale
%         rnoise = awgn(fadesig,SNR(n),'measured'); % Additive receiver noise, receiver snr
        rnoise = awgn(fadesig, SNR(n)); % additive noise with the assumption
                                        % that the received faded signal is
                                        % unit power (0dB).
%         rnoise = reshape(rnoise, [], nblock);
%         rnoise = rnoise(length_cp+1:end, :); % remove CP
        % for UW case, no need to remove anything
        %%
        % Channel estimation, ideal
        h = zeros(size(chan.pathGains,1), max(path_delays_norm)+1);
        h(:,path_delays_norm+1) = chan.pathGains;
        h_sq_norm = sum(abs(h(:)).^2); % channel pulse response energy for MFB calculation
        H = fft(h, nfft, 2);
        %%
        fdein = fft(rnoise,nfft); % convert into frequency domain for FDE
        ... FDE-LE to be inserted here
%         H_mmse = (H')./((abs(H).^2+ones(size(H))/(10^(SNR(n)/10))).'); %
%         mmse with receiver snr
        H_mmse = (H')./((abs(H).^2+ones(size(H))/(xtxPower/noisePower)).'); 
        % mmse with transmitted signal to receiver noise power ratio
        H_zf = (H')./((abs(H).^2).'); % zf
        fdeout_zf = fdein.*H_zf;
        fdeout_mmse = fdein.*H_mmse;
        xrx_zf = ifft(fdeout_zf,nfft); % convert back to time domain for detection
        xrx_mmse = ifft(fdeout_mmse,nfft); % convert back to time domain for detection
        % Demodulate in time domain
        z_zf=demodulate(modem.qamdemod(M_mod),xrx_zf);
        z_mmse=demodulate(modem.qamdemod(M_mod),xrx_mmse);
        %%
        % fd-dfe
%         kb = path_delays_norm(2:end); % non-zero indices correspond to the delays 
%                                       % (in symbol periods) of the feedback coefficients
        kb = (1: max(path_delays_norm)); % kb (full) non-zero indices
        V = zeros(length(kb), length(kb), nfft);
        v = zeros(length(kb), 1, nfft);
        for p = 1: nfft
            V(:,:,p) = (exp(1i*2*pi*kb*(p-1)/nfft)')*(exp(1i*2*pi*kb*(p-1)/nfft))/(abs(H(p)).^2+noisePower/xtxPower);
            v(:,:,p) = (exp(1i*2*pi*kb*(p-1)/nfft)')/(abs(H(p)).^2+noisePower/xtxPower);
%             for q = 1: length(kb)
%                 V(q,:,p) = exp(-1i*2*pi*kb(q)*(p-1)/nfft)*(exp(1i*2*pi*kb*(p-1)/nfft))/(abs(H(p)).^2+noisePower/xtxPower);
%                 v(q,:,p) = exp(-1i*2*pi*kb(q)*(p-1)/nfft)/(abs(H(p)).^2+noisePower/xtxPower);
%             end
        end
        V = sum(V,3); v = sum(v,3);
        V = squeeze(V); v = squeeze(v);
        fb = -V\v; % calculating non-zero feedback filter coefficients
        fbb = zeros(1, max(kb)); % concerning zero coefficients
        fbb(kb) = fb;
%         FB = 1 + fft(conj([0 fbb]),nfft); % convert to frequency domain
        FB = fft(conj([1 fbb]),nfft); % convert to frequency domain
        FB = FB.';
        ... feed foward filter to be inserted here
        W = ((H').*FB)./(noisePower/xtxPower+(abs(H.').^2));
        r = ifft(W.*fdein); % feedfoward filter output
        ... feed back structure
        xrx_dfe = zeros(size(xrx_zf));
        zrx_dfe = zeros(size(xrx_dfe));
        z_dfe = zeros(size(xrx_dfe));
        uw = xtx(end: -1: end-length(fbb)+1); % unique words for equalizer 
                                         % initialization, ideal cp adopted
                                         % here.
        for p = 1: length_pload
            zrx_dfe(p) = r(p) - conj(fbb)*uw;
            z_dfe(p) = demodulate(modem.qamdemod(M_mod),zrx_dfe(p)); % detection
%             xrx_dfe(p) = modulate(modem.qammod(M_mod),z_dfe(p)); % feedback symbol
%             xrx_dfe(p) = xrx_dfe(p)./sqrt(abs(xrx_dfe(p)).^2);
%             uw = [xrx_dfe(p); uw(1: end-1)]; % feedback filtering
            uw = [xtx(p); uw(1: end-1)]; % ideal feedback
        end
        %%
        ... IB-DFE to be inserted here
        n_ib = 4; % #of iterations
        rs = 0;
        z_ib = zeros(length_pload,1); % detected symbols
        x_ib = zeros(size(z_ib));
        Sd = zeros(nfft,1); % DFT of remodulated symbols
        xib_ext = zeros(size(Sd));
        C = zeros(nfft,1);  % DFT of feedforward filter
        B = zeros(nfft,1);  % DFT of feedback filter
        M_Sd = 1;  % Power of detected symbols, initial value set to non-zero to avoid conflits in denominator.
        for p = 1: n_ib
%             rs = Sd'*XTX/nfft;   % correlation of transmitted signal and
%                                  % detected signal. XTX need to be
%                                  % estimated at rx end in practice.
%             rs = xtx'*[x_ib;zeros(length_cp,1)]/(xtx'*xtx); % Chan & Wornell 2001 version
            ... forward filter update
            for q = 1: nfft
%                 C(q) = H(q)'/(nfft*noisePower+abs(H(q)).^2*(1-abs(rs).^2/(M_XTX*M_Sd))*M_XTX);
                C(q) = H(q)'/(nfft*noisePower+xtx'*xtx*(1-abs(rs).^2)*(abs(H(q)).^2)); % Chan & Wornell 2001 version
            end
            gamma = H*C/nfft;
%             C = M_XTX*(1-abs(rs).^2*gamma/(M_XTX*M_Sd))*C;            
            ... feedback filter update
            for q = 1: nfft
%                 B(q) = -(rs*(H(q)*C(q)-gamma))/M_Sd;
                B(q) = -rs*(H(q)*C(q)-gamma); % Chan & Wornell 2001 version
            end
            U = fdein.*C+Sd.*B; % combine with feedback symbols
            u = ifft(U,nfft);
            zrx_ib = u(1:length_pload);
            z_ib = demodulate(modem.qamdemod(M_mod),zrx_ib); % hard detection in time domain
%             % tp
%             if (rs < 1 && rs > 0)
%                 symerr(z_mmse(1:length_pload),x)
%                 symerr(z_dfe(1:length_pload),x)
%                 symerr(z_ib(1:length_pload),x)
%             end
%             % end of tp
            x_ib = modulate(modem.qammod(M_mod),z_ib);
            xib_ext = [x_ib;zeros(length_cp,1)];
            xib_ext = xib_ext/sqrt(sum(abs(xib_ext).^2)/nfft);
            Sd = fft([x_ib;zeros(length_cp,1)],nfft);
%             Sd = XTX; % ideal feedback, for test purpose only
            M_Sd = Sd'*Sd/nfft;  % Power of detected symbols
            rs = xtx'*xib_ext/(sqrt(xtx'*xtx)*sqrt(xib_ext'*xib_ext)); % Chan & Wornell 2001 version
%             rs = xtx'*xib_ext/(xtx'*xtx);
        end
        %%
        % Performance evaluation
        [numzf]  = symerr(z_zf(1:length_pload),x);   % symbol error rate
        [nummse] = symerr(z_mmse(1:length_pload),x); % symbol error rate
        [numdfe] = symerr(z_dfe(1:length_pload),x);  % symbol error rate
        [numib]  = symerr(z_ib(1:length_pload),x);   % symbol error rate
        % Compute MFB
%         snr_mfb = (10^(SNR(n)/10))*h_sq_norm; % with receiver snr
        snr_mfb = (xtxPower*h_sq_norm/noisePower); % with transmitted signal to receiver noise
                                         % power ratio
        kap = 1; % kap = (.5*d_min)^2/((i^2+q^2)/N); where d_min is the 
                 % minimum distance between a constellation pair, and N
                 % equals to 2 when using quadrature (IQ) modulation.
        mfb = snr_mfb*kap;
        ber_mfb(n) = ber_mfb(n) + .5*erfc(sqrt(mfb)/sqrt(2)); % Q(x) = .5*erfc(x/sqrt(2))
        numerr(:,n) = numerr(:,n)+[numzf;nummse;numdfe;numib];
        waitbar(((n-1)*numTrials+m)/(numTrials*length(SNR)))
    end
end
close(h_waitbar)
BER = numerr/(nfft*nblock*numTrials);
ber_mfb = ber_mfb/numTrials;

%%
% Plot BER results.
figure();
semilogy(SNR,BER(1,:),'-*',SNR,BER(2,:),'-o',SNR,BER(3,:),'-p',SNR,BER(4,:),'-s',SNR,ber_mfb,'r-');
legend('ZF','MMSE','DFE','IB','MFB');
xlabel('SNR (dB)'); ylabel('BER');
title('QPSK-SC-FDE over Rayleigh Fading Channel');
ylim([5e-6 1]);
grid on
%%

filename = 'sc_fde_ib_dfe_chan_wornell.mat';
save(filename, 'SNR', 'BER', 'ber_mfb', 'path_delays_norm', 'doppler_max', 'pdp');

% End of script
% [1] A. Chan and G. Wornell, ¡°A class of block-iterative equalizers for intersymbol interference channels: fixed channel results,¡± Communications, IEEE Transactions on,  vol. 49, 2001, pp. 1966-1976.
% [2] N. Benvenuto and S. Tomasin, ¡°Iterative design and detection of a DFE in the frequency domain,¡± Communications, IEEE Transactions on,  vol. 53, 2005, pp. 1867-1875.
