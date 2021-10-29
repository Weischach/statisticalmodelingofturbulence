function p1 = psdperturbseries(b, f1, r1, tau, variance, coeffs)
rtemp = 1;
R = zeros(size(r1));
for k = 1:coeffs
    rtemp = rtemp.*r1;
    R = R + b(k)*rtemp;
end
for kk = 1:length(f1)
    pwwtemp = R.*cos(2*pi*f1(kk).*tau);
    p1(kk) = 4*trapz(tau, pwwtemp)*variance;
end
end