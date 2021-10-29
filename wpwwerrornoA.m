function error = wpwwerrornoA(Betas, p1, f1, al_coeff, T, t, variance, coeffs)
alphas = dot(Betas,al_coeff);
tuu = (alphas/T)*t;
r1 = ((2^(2/3))/gamma(1/3))*[(tuu.^(1/3)).*besselk(1/3,tuu) - 0.5*(tuu.^(4/3)).*besselk(2/3,tuu)];
r1(1) = 1;
rtemp = 1;
Ruu = zeros(size(r1));
for k = 1:coeffs
    rtemp = rtemp.*r1;
    Ruu = Ruu + Betas(k)*rtemp;
end
for kk = 1:length(f1)
    pwwtemp = Ruu.*cos(2*pi*f1(kk).*t);
    p2(kk) = 4*trapz(t, pwwtemp)*variance;
end
p2 = f1.*reshape(p2,size(f1));
p1 = f1.*reshape(p1,size(f1));
s = length(p2);
kt = min(p2);
kt2 = (kt^2) - (abs(kt)*kt);
for k = 1:s
    ertemp(k) = ((p2(k)-p1(k)))^2;
end
error = sum(ertemp) + (10^20)*kt2;
% error = sum(ertemp);
% fprintf('The error is %f \n', error)
end