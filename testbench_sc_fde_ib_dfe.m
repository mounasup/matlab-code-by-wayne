% File: testbench_sc_fde_ib_dfe.m
% -------------------------------
% This script This script is a testbench containning debug function calls 
% to the sc_fde_ib_dfe.m

close all
clear all

% dbstop in sc_fde_ib_dfe at 157 if (rs < XTX'*XTX/nfft && rs > 0)
dbstop in sc_fde_ib_dfe at 186 if (rs < 1 && rs > 0)
dbstop in sc_fde_ib_dfe at 203 if (rs < 1 && rs > 0)

sc_fde_ib_dfe

% End of script