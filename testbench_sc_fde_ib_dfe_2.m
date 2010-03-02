% File: testbench_sc_fde_ib_dfe_2.m
% ---------------------------------
% This script tests IB-DFE implementation [1] with a specified work space.

close all
clear all
clc
% load testbench_sc_fde_ib_dfe_.9669-.0331i.mat
% load testbench_sc_fde_ib_dfe_.9587-.0083i.mat
load testbench_sc_fde_ib_dfe_.8678-.0165i.mat
% load testbench_sc_fde_ib_dfe_.9091-.0083i.mat

% IB-DFE process
        n_ib = 10; % #of iterations
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
                                 % detected signal. XTX need to be
                                 % estimated at rx end in practice.
%             rs = xtx'*[x_ib;zeros(length_cp,1)]/(xtx'*xtx); % Chan & Wornell 2001 version
            ... forward filter update
            for q = 1: nfft
%                 C(q) = H(q)'/(nfft*noisePower+abs(H(q)).^2*(1-abs(rs).^2/(M_XTX*M_Sd))*M_XTX);
%                 C(q) = H(q)'/(nfft*noisePower+xtx'*xtx*(1-abs(rs).^2)*(abs(H(q)).^2)); % Chan & Wornell 2001 version
                C(q) = (xtx'*xtx)*H(q)'/(nfft*noisePower+(xtx'*xtx)*(1-abs(rs).^2)*(abs(H(q)).^2)); % Chan & Wornell 2001 version
            end
            gamma = H*C/nfft;
%             gamma = H*C;
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
            % tp
            if (rs < 1 && rs > 0)
                symerr(z_mmse(1:length_pload),x)
                symerr(z_dfe(1:length_pload),x)
                symerr(z_ib(1:length_pload),x)
            end
            % end of tp
            x_ib = modulate(modem.qammod(M_mod),z_ib);
            xib_ext = [x_ib;zeros(length_cp,1)];
            xib_ext = xib_ext/sqrt(sum(abs(xib_ext).^2)/nfft);         
            Sd = fft(xib_ext,nfft);
%             Sd = XTX; % ideal feedback, for test purpose only
%             M_Sd = Sd'*Sd/nfft;  % Power of detected symbols
            rs = xtx'*xib_ext/(sqrt(xtx'*xtx)*sqrt(xib_ext'*xib_ext))
%             rs = xtx'*xib_ext/(xtx'*xtx)
        end

% End of script
% [1] A. Chan and G. Wornell, ¡°A class of block-iterative equalizers for
% intersymbol interference channels: fixed channel results,¡±
% Communications, IEEE Transactions on,  vol. 49, 2001, pp. 1966-1976.