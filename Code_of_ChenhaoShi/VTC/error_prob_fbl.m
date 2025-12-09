function err = error_prob_fbl(snr, m, r)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

V=1-1./(snr+1).^2;
w=(m./V).^0.5.*(log2(1+snr)-r);

% disp(imag(w)~=0);
w(imag(w)~=0)=-inf;
w(snr<0)=-inf;
%w=real(w)+log(double(imag(w)==0))+log(double(snr>=0));

err=normcdf(-w);

end