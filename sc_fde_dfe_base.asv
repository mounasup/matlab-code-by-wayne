% File: sc_fde_dfe_base.m
% -----------------------
% This script implements SC-FDE-DFE system based on 'sc_fde_le_base.m'
% and [1][2]. It uses receiver snr to for the evaluation of ber and 
% ber_mfb perfomance.
close all
clear all

%%
BW = 2e6; % bandwidth in bits/sec.
Ts = 1/BW; % symbol interval
path_delays_norm = [0 2 5 7]; % multipath delay spread normalized to Ts
doppler_max = 0; % max doppler shift in Hz
pdp = [0 0 -5 -5]; % power delay profile

nfft = 64; % #of symbols in a fde block
length_cp = max(path_delays_norm); % #of symbols in cyclic prefix
nblock = 1; % #of fde blocks in each frame, >1 consider IBI; <=1 without considering IBI
numTrials = 100; % #of frames in each evaluation
M_mod = 4; % Alphabet size
%%
% Compute error rate for different values of SNR.
SNR = 0:2:20; % Range of SNR values, in dB. receiver snr
SNR = 16;
numerr = zeros(3,length(SNR)); % zf, mmse, dfe
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
        x = randint(nfft,nblock,M_mod);
        xtx=modulate(modem.qammod(M_mod),x);
        xtxPower = sum(abs(xtx).^2,1)/size(xtx,1); % check transmitted signal power along each tx antenna
        xtx = xtx*diag(1./sqrt(xtxPower))./sqrt(1); % normalized and fixed total power
        xtx = [xtx(end-length_cp+1: end, :); xtx]; % insert CP
        xtx = reshape(xtx, [], 1);
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
        rnoise = reshape(rnoise, [], nblock);
        rnoise = rnoise(length_cp+1:end, :); % remove CP
        %%
        % Channel estimation, ideal
        h = zeros(size(chan.pathGains,1), max(path_delays_norm)+1);
        h(:,path_delays_norm+1) = chan.pathGains;
        h_sq_norm = sum(abs(h(:)).^2); % channel pulse response energy for MFB calculation
        H = fft(h, nfft, 2);
        %%
        fdein = fft(rnoise); % convert into frequency domain for FDE
        ... FDE-LE to be inserted here
%         H_mmse = (H')./((abs(H).^2+ones(size(H))/(10^(SNR(n)/10))).'); %
%         mmse with receiver snr
        H_mmse = (H')./((abs(H).^2+ones(size(H))/(xtxPower/noisePower)).'); 
        % mmse with transmitted signal to receiver noise power ratio
        H_zf = (H')./((abs(H).^2).'); % zf
        fdeout_zf = fdein.*H_zf;
        fdeout_mmse = fdein.*H_mmse;
        xrx_zf = ifft(fdeout_zf); % convert back to time domain for detection
        xrx_mmse = ifft(fdeout_mmse); % convert back to time domain for detection
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
        for p = 1: nfft
            zrx_dfe(p) = r(p) - conj(fbb)*uw;
            z_dfe(p) = demodulate(modem.qamdemod(M_mod),zrx_dfe(p)); % detection
            xrx_dfe(p) = modulate(modem.qammod(M_mod),z_dfe(p)); % feedback symbol
            xrx_dfe(p) = xrx_dfe(p)./sqrt(abs(xrx_dfe(p)).^2);
            uw = [xrx_dfe(p); uw(1: end-1)]; % feedback filtering
%             uw = [xtx(length_cp+p); uw(1: end-1)]; % ideal feedback
        end
        %%
        % Performance evaluation
        [numzf] = symerr(z_zf,x); % symbol error rate
        [nummse] = symerr(z_mmse,x); % symbol error rate
        [numdfe] = symerr(z_dfe,x); % symbol error rate
        % Compute MFB
%         snr_mfb = (10^(SNR(n)/10))*h_sq_norm; % with receiver snr
        snr_mfb = (xtxPower*h_sq_norm/noisePower); % with transmitted signal to receiver noise
                                         % power ratio
        kap = 1; % kap = (.5*d_min)^2/((i^2+q^2)/N); where d_min is the 
                 % minimum distance between a constellation pair, and N
                 % equals to 2 when using quadrature (IQ) modulation.
        mfb = snr_mfb*kap;
        ber_mfb(n) = ber_mfb(n) + .5*erfc(sqrt(mfb)/sqrt(2)); % Q(x) = .5*erfc(x/sqrt(2))
        numerr(:,n) = numerr(:,n)+[numzf;nummse;numdfe];
        waitbar(((n-1)*numTrials+m)/(numTrials*length(SNR)))
    end
end
close(h_waitbar)
BER = numerr/(nfft*nblock*numTrials);
ber_mfb = ber_mfb/numTrials;

%%
% Plot BER results.
figure();
semilogy(SNR,BER(1,:),'-*',SNR,BER(2,:),'-o',SNR,BER(3,:),'-p',SNR,ber_mfb,'r-');
legend('ZF','MMSE','DFE','MFB');
xlabel('SNR (dB)'); ylabel('BER');
title('QPSK-SC-FDE over Rayleigh Fading Channel');
ylim([5e-6 1]);
grid on
%%

filename = 'sc_fde_le_base.mat';
save(filename, 'SNR', 'BER', 'ber_mfb', 'path_delays_norm', 'doppler_max', 'pdp');

% End of script
% [1] D. Falconer, S. Ariyavisitakul, A. Benyamin-Seeyar, and B. Eidson, ¡°Frequency domain equalization for single-carrier broadband wireless systems,¡± Communications Magazine, IEEE,  vol. 40, 2002, pp. 58-66.
% [2] N. Benvenuto and S. Tomasin, ¡°On the comparison between OFDM and single carrier modulation with a DFE using a frequency-domain feedforward filter,¡± Communications, IEEE Transactions on,  vol. 50, 2002, pp. 947-955.

