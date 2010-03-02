% File: mimochan_filter.m
% -----------------------
% This function implements mimo channel filtering, using filter function for
% each individual link.

function [fadesig h h_sq_norm] = mimochan_filter(chan_mimo, xtx)
[nr nt] = size(chan_mimo);
path_delays_norm = chan_mimo{1,1}.PathDelays/chan_mimo{1,1}.InputSamplePeriod;
fadesig = zeros(nr,nt,size(xtx,1)); % fading signal for each link
h = zeros(nr,nt,max(path_delays_norm)+1); % channel impulse response for each link
for m = 1: nr
    for n = 1: nt
        fadesig(m,n,:) = filter(chan_mimo{m,n},xtx(:,n));
        h(m,n,path_delays_norm+1) = chan_mimo{m,n}.PathGains(1,:);
    end
end
fadesig = squeeze(sum(fadesig,2)).';
fadesig = reshape(fadesig,[],nr);
h_sq_norm = sum(abs(h).^2,3); % channel pulse response energy

% End of function