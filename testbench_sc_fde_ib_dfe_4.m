% File: testbench_sc_fde_ib_dfe_4.m
% ---------------------------------
% This script test IB-DFE performance with respect to various block length, 
% fixed channel length, and fixed SNR, based on sc_fde_ib_dfe_func.m

close all
clear all
clc

pd = (0:1:31);
nfft = 2.^(5:12);
snr = 12;
ntrial = 100000;
ber = zeros(4,length(nfft));
mfb = zeros(length(nfft));

for n = 1: length(nfft)
    nfft(n)
    [ber(:,n) mfb(n)] = sc_fde_ib_dfe_func(pd, nfft(n), snr, ntrial);
end

% ber
% mfb

filename = 'testbench_sc_fde_ib_dfe_4.mat';
save(filename, 'pd', 'nfft', 'snr', 'ber', 'mfb');
% End of script