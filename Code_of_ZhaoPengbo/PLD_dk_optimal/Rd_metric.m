function rd = Rd_metric(eBM, eBK, eEM, eEK)
% 有效欺骗率 R_d 指标
% 与 OFDM 论文一致：
%   R_d = [1 - (1-eBM)*eBK] * (1-eEM) * eEK

rd = (1 - (1 - eBM) .* eBK) .* (1 - eEM) .* eEK;
end