% File: sc_fde_dfe_mimo_base.m
% ----------------------------
% This script implements SC-FDE-DFE-MIMO system based on 'sc_fde_dfe_base.m'
% and [1][2].
close all
clear all

%%
BW = 2e6; % bandwidth in bits/sec.
Ts = 1/BW; % symbol interval
path_delays_norm = [0 2 5 7]; % multipath delay spread normalized to Ts
doppler_max = 0; % max doppler shift in Hz
                 % the bits that are not affected by doppler can be computed
                 % like this:
                 % bit_rate/(100*fd), where fd = c/(v*fc)
                 % e.g. 2e6/(100*1) = 2e4
                 % thus 20000 consecutive bits will not be affected by
                 % doppler, and choose no more than 20000 within each tial
pdp = [0 0 -5 -5]; % power delay profile, in dB
nt = 4; % #of tx antennas
nr = 4; % #of rx antennas

nfft = 64; % #of symbols in a fde block
length_cp = max(path_delays_norm); % #of symbols in cyclic prefix
length_pload = nfft-length_cp; % #of play load symbols in null UW scheme
numTrials = 10000; % #of frames in each evaluation
M_mod = 4; % Alphabet size
%%
% Compute error rate for different values of SNR.
SNR = 8:2:26; % Range of SNR values, in dB. receiver snr
numerr = zeros(3,nt,length(SNR)); %       column-wise ... 
                                    % zf
                                    % mmse 
                                    % dfe
ber_mfb = zeros(nt,length(SNR));
h_waitbar = waitbar(0,'Please wait...');
for n = 1:length(SNR)
    % mimo channel initialization, refers ot mimochan_init.m for details
    chan_mimo = mimochan_init(nt,nr,Ts,doppler_max,path_delays_norm,pdp);
%     chan = mimochan(nt,nr,Ts,doppler_max,path_delays_norm*Ts,pdp); % mimo channel method
% %     chan = rayleighchan(Ts,doppler_max,path_delays_norm*Ts,pdp);
%     chan.StorePathGains = 1; % the complex path gain vector is stored 
%                                 % as the channel filter function processes 
%                                 % the signal
%     chan.ResetBeforeFiltering = 1; % each call to filter resets the 
%                                       % state of chan before filtering
%     chan.NormalizePathGains = 1;  % If 1, the Rayleigh fading process 
%                                       % is normalized such that the expected 
%                                       % value of the path gains' total power is 1.
%%
    for m = 1:numTrials % quasi-static channel modeling
        % Create a random digital message
        % Use M-QAM modulation.
%         x = randi([0,M_mod-1], nfft, nt); % for CP scheme
        x = randi([0,M_mod-1], length_pload, nt); % for UW scheme
%         x = randint(length_pload,nt,M_mod); % for old version MATLAB
%         compatible
        xtx_pload = modulate(modem.qammod(M_mod),x);
        xtxPower = sum(abs(xtx_pload).^2,1)/size(xtx_pload,1); % check transmitted signal power along each tx antenna
        xtx_pload = xtx_pload*diag(1./sqrt(xtxPower))./sqrt(nt); % normalized and fixed total power
%         xtx = [xtx_pload(end-length_cp+1: end, :); xtx_pload]; % insert CP
        xtx = [xtx_pload; zeros(length_cp,nt)]; % add UW at the end of a block.
                                          % the UW at the beginning of a
                                          % block is omitted, since the UW
                                          % we choosed here is a null
                                          % sequence, it will not interfere
                                          % with the pay load.
        % mimo channel filtering, refers to mimochan_filter.m for details
        [fadesig h h_sq_norm] = mimochan_filter(chan_mimo, xtx);
%         fadesig = filter(chan, xtx); % Effect of channel, quasi-static
%         % Caution! The output is incorrect. Problem to be fixed here:
%         % When max_doppler is set to 0, there will be an incorrect
%         % output from the filter(mimo) function. This may be a bug of the
%         % mimochan module. When max_doppler is set > 0, the output seems to
%         % be OK, but the final performance (MIMO 1x1 scenario) is not as 
%         % good as the SISO scenario from rayleighchan module. It is
%         % still not clear that whether or not other causes (e.g. incorrect 
%         % snr interpretation) exist.
%         h = zeros(nr,nt,max(path_delays_norm)+1);
%         for p = 1: nr
%             for q = 1: nt
%                 h(p,q,path_delays_norm+1) = chan.PathGains(1,:,q,p); % assume quasi-static channel response, although doppler is not zero
%             end
%         end
%         h_sq_norm = sum(abs(h).^2,3); % channel pulse response energy
        %%
        % energy/power re-calculation inserted here...
        xtxPower = sum(abs(xtx).^2,1)/size(xtx,1); % check transmitted signal power along each tx antenna
        xtxPower_dB = 10*log10(xtxPower); % convert to dB
        xtxPower_avg = sum(xtxPower)/size(xtx,2); % averaged tx power, with path gain (power) normalized to 1
        xtxPower_dB_avg = 10*log10(xtxPower_avg); % convert to dB      
%         fadesigPower_mf = h_sq_norm*(xtxPower.'); % matched filter output signal at each rx antenna
%         fadesigPower_mf_dB = 10*log10(fadesigPower_mf);
%         noisePower_mf_dB = fadesigPower_mf_dB - SNR(n); % noise power at matched filter output
%         noisePower_mf = 10.^(noisePower_mf_dB/10);        
        fadesigPower = sum(abs(fadesig).^2,1)/size(fadesig,1); % check fading signal power
        fadesigPower_dB = 10*log10(fadesigPower); % convert to dB
        % spatial noise whitening
%         noisePower_dB = fadesigPower_dB-SNR(n); % compute noise power with respect to 
%                                                 % the instanteneous faded signal power, 
%                                                 % in dB
        noisePower_dB = zeros(size(fadesigPower_dB))-SNR(n); % compute noise power, with 
                                                  % the assumption that the
                                                  % received signal at each
                                                  % rx-antenna is unit
                                                  % power (0dB).
        noisePower = 10.^(noisePower_dB/10); % convert back to linear scale
        snow = diag(sqrt(max(noisePower)./noisePower));
        noisePower_snow = max(noisePower);
        noisePower_snow_dB = 10*log10(noisePower_snow); % convert to dB
        % awgn channel
        rnoise = zeros(size(fadesig));
        for p = 1: nr
%             rnoise(:,p) = awgn(fadesig(:,p), SNR(n), fadesigPower_dB(p));
%             rnoise(:,p) = awgn(fadesig(:,p), SNR(n), 'measured');
            rnoise(:,p) = awgn(fadesig(:,p), SNR(n)); % additive noise with the assumption
                                                      % that the received
                                                      % faded signal is
                                                      % unit power (0dB).
        end
        rnoise = rnoise * snow; % noise whitening
%         rnoise = rnoise(length_cp+1:end, :); % remove CP
        % for UW case, no need to remove anything
        %%
        % Channel estimation, ideal
        H = fft(h, nfft, 3);
        %%
        fdein = fft(rnoise,nfft,1).'; % convert into frequency domain for FDE
        ... MIMO FDE-LE to be inserted here
        W_mmse = zeros(nt,nr,nfft);
        R_mmse = zeros(nt,nfft);
        for p = 1: nfft
            W_mmse(:,:,p) = (H(:,:,p)')/(H(:,:,p)*(H(:,:,p)')+... % mmse
                (noisePower_snow./xtxPower_avg)*eye(nr));
            R_mmse(:,p) = W_mmse(:,:,p)*fdein(:,p);
        end
        W_zf = zeros(nt,nr,nfft);
        R_zf = zeros(nt,nfft);
        for p = 1: nfft
            W_zf(:,:,p) = pinv(H(:,:,p)); % zf
            R_zf(:,p) = W_zf(:,:,p)*fdein(:,p);
        end
        % Demodulate in time domain
        r_mmse = ifft(R_mmse, nfft, 2); % convert back to time domain for detection
        r_zf = ifft(R_zf, nfft, 2); % convert back to time domain for detection
        z_zf=demodulate(modem.qamdemod(M_mod),r_zf.');
        z_mmse=demodulate(modem.qamdemod(M_mod),r_mmse.');
        %%
        ... MIMO FD-DFE to be inserted here
        kb = (1: max(path_delays_norm)); % kb non-zero indices
        V = zeros(nt*length(kb), nt*length(kb), nfft);
        v = zeros(nt*length(kb), nt*1, nfft);
        for p = 1: nfft
            V(:,:,p) = kron((exp(1i*2*pi*kb*(p-1)/nfft)')*(exp(1i*2*pi*kb*(p-1)/nfft)),...
                pinv((H(:,:,p)'*H(:,:,p)+(noisePower_snow./xtxPower_avg)*eye(nt)))); % for uncorrelated signal and noise
            v(:,:,p) = kron((exp(1i*2*pi*kb*(p-1)/nfft)'),...
                pinv((H(:,:,p)'*H(:,:,p)+(noisePower_snow./xtxPower_avg)*eye(nt))));
        end
        V = sum(V,3); v = sum(v,3);
%         q_core = zeros(nt,nt,nfft); % another way for computing V & v
%         for p = 1: nfft
%             q_core(:,:,p) = inv(H(:,:,p)'*H(:,:,p)+(noisePower_snow./xtxPower_avg)*eye(nt));
%         end
%         Q = ifft(q_core,nfft,3);
%         V = zeros(nt*length(kb), nt*length(kb));
%         v = zeros(nt*length(kb), nt*1);
%         Q_kb = Q(:,:,1:max(kb)); % corresponding to index 0: kb-1
%         for p = kb
%             V((p-1)*nt+1:p*nt,:) = reshape(Q_kb,nt,[]);
%             Q_kb = cat(3,Q(:,:,p+1)',Q_kb);
%             Q_kb = Q_kb(:,:,1:max(kb));
%             v((p-1)*nt+1:p*nt,:) = Q(:,:,p+1)'; % corresponding to index 1: kb
%         end                
        fb = -V\v; % calculating non-zero feedback filter coefficients
        ... fb to fbb to FB convertion refinement inserted here
        fbb = zeros(nt,nt,max(kb)); % concerning zero coefficients
        fbb_conj = zeros(size(fbb));
        for p = kb
            fbb(:,:,p) = fb((p-1)*nt+1:p*nt,:);
            fbb_conj(:,:,p) = fbb(:,:,p)'; % fbb is not hermitian, so complex conjugate needed
        end
        FB = fft(cat(3,eye(nt), fbb_conj),nfft,3); % convert to frequency domain
        ... feed foward filter to be inserted here
        W_dfe = zeros(nt,nr,nfft);
        R = zeros(nt,nfft);
        for p = 1: nfft
            W_dfe(:,:,p) = FB(:,:,p)*(H(:,:,p)')/(H(:,:,p)*(H(:,:,p)')+...
                (noisePower_snow./xtxPower_avg)*eye(nr));
            R(:,p) = W_dfe(:,:,p)*fdein(:,p);
        end
        r = ifft(R,nfft,2); % feedfoward filter output     
        ... feed back structure
        xrx_dfe = zeros(size(r));
        zrx_dfe = zeros(size(xrx_dfe));
        z_dfe = zeros(size(xrx_dfe));
%         uw = xtx(end: -1: end-size(fbb,3)+1,:).'; % unique words for equalizer 
%                                          % initialization, ideal cp adopted
%                                          % here.
        uw = zeros(nt,size(fbb,3)); % for null uw case
        for p = 1: length_pload
            zrx_dfe(:,p) = r(:,p) - reshape(fbb_conj,nt,[])*reshape(uw,[],1);
            z_dfe(:,p) = demodulate(modem.qamdemod(M_mod),zrx_dfe(:,p).').'; % detection
%             xrx_dfe(:,p) = modulate(modem.qammod(M_mod),z_dfe(:,p).').'; % feedback symbol
%             % xrx_dfe should be power normalized
%             xrx_dfe(:,p) = xrx_dfe(:,p)./(sqrt(abs(xrx_dfe(:,p)).^2)*sqrt(nt));
%             uw = [xrx_dfe(:,p) uw(:,1: end-1)]; % feedback filtering
%             uw = [xtx(length_cp+p,:).' uw(:,1: end-1)]; % ideal feedback for non-zero cp case
            uw = [xtx(p,:).' uw(:,1:end-1)]; % ideal feedback for null uw case
        end
        %%
        % Performance evaluation
        [numzf] = symerr(z_zf(1:length_pload,:),x,'column-wise'); % symbol error rate
        [nummse] = symerr(z_mmse(1:length_pload,:),x,'column-wise'); % symbol error rate
        [numdfe] = symerr(z_dfe(:,1:length_pload).',x,'column-wise'); % symbol error rate
        numerr(:,:,n) = numerr(:,:,n)+[numzf;nummse;numdfe];
        %%
        % Compute MFB
        % need for refining...
        xtxPower_mfb = (xtxPower).*sum((snow.^2)*h_sq_norm,1); % MF signal power for MFB
        xtxPower_mfb_dB = 10*log10(xtxPower_mfb);
        snr_mfb = (xtxPower_mfb/noisePower_snow); % with transmitted signal to receiver noise
                                                  % power ratio
        snr_mfb_dB = 10*log10(snr_mfb);
        kap = 1; % kap = (.5*d_min)^2/((i^2+q^2)/N); where d_min is the 
                 % minimum distance between a constellation pair, and N
                 % equals to 2 when using quadrature (IQ) modulation.
        mfb = snr_mfb*kap;
        ber_mfb(:,n) = ber_mfb(:,n) + .5*erfc(sqrt(mfb.')/sqrt(2)); % Q(x) = .5*erfc(x/sqrt(2))
        %%
        waitbar(((n-1)*numTrials+m)/(numTrials*length(SNR)))
    end
end
close(h_waitbar)
BER = numerr/((nfft-length_cp)*numTrials);
BER_oval = squeeze(sum(BER,2)/size(BER,2));
ber_mfb = ber_mfb/numTrials;
ber_mfb_oval = squeeze(sum(ber_mfb,1)/size(ber_mfb,1));

%%
% Plot BER results.
figure();
semilogy(SNR,BER_oval(1,:),'-*',SNR,BER_oval(2,:),'-o',SNR,BER_oval(3,:),'-p',SNR,ber_mfb_oval,'r-');
legend('ZF','MMSE','DFE','MFB');
xlabel('SNR (dB)'); ylabel('BER');
title('QPSK-SC-FDE over Rayleigh Fading MIMO Channel');
ylim([1e-6 max(BER_oval(:))]);
grid on
%%

filename = 'sc_fde_dfe_mimo_base.mat';
save(filename, 'nt', 'nr', 'SNR', 'BER', 'BER_oval', 'ber_mfb', 'ber_mfb_oval', 'path_delays_norm', 'doppler_max', 'pdp');

% End of script
% [1] J. Tubbax, L. Van der Perre, S. Donnay, and M. Engels, ¡°Single-carrier communication using decision-feedback equalization for multiple antennas,¡± Communications, 2003. ICC '03. IEEE International Conference on, 2003, pp. 2321-2325 vol.4.
% [2] A. Feng, Q. Yin, and J. Fan, ¡°Hybrid Two-Stage Decision-Feedback
% Equalization for Single-Carrier Multiple-Input Multiple-Output Systems,¡± IEICE TRANS. COMMUN.,  vol. E92-B, Jul. 2009.