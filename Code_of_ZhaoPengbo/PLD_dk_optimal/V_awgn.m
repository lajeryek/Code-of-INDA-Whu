function v = V_awgn(gamma)
% AWGN 信道色散 V(SNR)
v = 1 - 1 ./ (1 + gamma).^2;
end