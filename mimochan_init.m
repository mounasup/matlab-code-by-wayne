% File: mimochan_init.m
% ---------------------
% This function initialize the MIMO channel, using rayleighchan for each
% individual link.
function [chan_mimo] = mimochan_init(nt,nr,Ts,doppler_max,path_delays_norm,pdp)
chan_mimo = cell(nr,nt); % creating cell arrays for rayleigh channel objects.
for m = 1: nr
    for n = 1: nt
        chan_mimo{m,n} = rayleighchan(Ts,doppler_max,path_delays_norm*Ts,pdp);
        chan_mimo{m,n}.StorePathGains = 1; % the complex path gain vector is stored 
                                 % as the channel filter function processes 
                                 % the signal
        chan_mimo{m,n}.ResetBeforeFiltering = 1; % each call to filter resets the 
                                       % state of chan before filtering
        chan_mimo{m,n}.NormalizePathGains = 1;   % If 1, the Rayleigh fading process 
                                       % is normalized such that the expected 
                                       % value of the path gains' total
                                       % power is 1.
    end
end

% End of function
