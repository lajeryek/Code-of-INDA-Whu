function [e] = error(m_n,SNR_normal,d)
   V_r = 1-(1+SNR_normal).^(-2);%信道色散 
   C_r=log2(1+SNR_normal);%香农容量
   e = qfunc(((m_n/V_r).^(1/2)).*(C_r-(d./m_n)));%erro_Q函数
   if m_n==0
       e=1;
   end
end
