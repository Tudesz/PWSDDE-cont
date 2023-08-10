function plot_spectrum(mu)
%PLOT_SPECTRUM Plot characteristic multipliers on Re-Im plane
% Input:
%   mu: characteristic multipliters

scatter(real(mu),imag(mu),50,'gx'); 
hold on
scatter(real(mu(abs(mu)>1)),imag(mu(abs(mu)>1)),50,'rx'); 
plot(real(exp(1i*(0:pi/100:2*pi))), imag(exp(1i*(0:pi/100:2*pi))),':k');
hold off; pbaspect([1 1 1]); axis('equal'); box on;
xlabel('$\Re$'); ylabel('$\Im$');
end

