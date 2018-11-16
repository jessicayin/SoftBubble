function y = generate_gaussian_samples(mu, sigma)
% Use The Box–Muller transform.
u1 = rand();
u2 = rand();

z1 = sqrt(-2*log(u1)) * cos(2*pi*u2);
%z2 = sqrt(-2*log(u1)) * sin(2*pi*u2);

y = mu + sigma * z1; % Drop z2.

end