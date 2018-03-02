clear all;

floor_var = 0;
% Define two Gaussians
p1 = 0.3; p2 = 1-p1;
mu1 = -0.8; mu2 = 1.2;
s1 = 0.52; s2 = 0.35;

% Draw samples
N = 50;
samples = [randn(N*p1,1)*s1+mu1;randn(N*p2,1)*s2+mu2];

px = exp(- ([-6:0.1:6] - mu1).^2/2/s1/s1)/sqrt(2*pi)/s1 * p1 ...
+ exp(- ([-6:0.1:6] - mu2).^2/2/s2/s2)/sqrt(2*pi)/s2 * p2;

[a, ax] = hist(samples, 10);
a = a/N;
clf;
bar(ax, a);
set(gca, 'Xlim', [-6 6], 'FontSize', 18);
hold;
plot([-6:0.1:6], px, 'r')

% starting EM, giving initial values
ep1  = 0.2; ep2 = 1-ep1;
emu1 = -2; emu2 = 6;
es1 =1; es2 = 1;

samples = sort(samples);
for i = 1:200
    
    pxx = exp(- ([-6:0.1:6] - emu1).^2/2/es1/es1)/sqrt(2*pi)/es1 * ep1 ...
        + exp(- ([-6:0.1:6] - emu2).^2/2/es2/es2)/sqrt(2*pi)/es2 * ep2;

    clf;
    bar(ax, a);
    set(gca, 'Xlim', [-6 6], 'FontSize', 18);
    hold on;
    plot([-6:0.1:6], px, 'r')
    plot([-6:0.1:6], pxx, 'b')
    pause
    % posterior

    % computing in log
    logpx1 = - (samples - emu1).^2/2/es1/es1 - log (sqrt(2*pi)*es1) + log(ep1);
    logpx2 =  - (samples - emu2).^2/2/es2/es2 - log(sqrt(2*pi)*es2) + log(ep2);
    maxx = max(logpx1, logpx2);
    logpx1 = logpx1 - maxx;
    logpx2 = logpx2 - maxx;
    like = exp(logpx1)+exp(logpx2);
    post1 = exp(logpx1)./like;
    fprintf('ll=%.2e ', sum(log(like)+maxx));


    % update
    N1 = sum(post1); N2 = sum(1-post1);
    ep1 = N1/(N1+N2); ep2 = 1-ep1;
    emu1 = sum( post1 .* samples)/N1;
    emu2 = sum( (1-post1) .* samples)/N2;
    if floor_var ==1 
        es1 = sqrt(sum( post1.* (samples-emu1).^2)/N1)+0.1*s1;
        es2 = sqrt(sum( (1-post1).* (samples-emu2).^2)/N2)+0.1*s2;
    else
        es1 = sqrt(sum( post1.* (samples-emu1).^2)/N1);
        es2 = sqrt(sum( (1-post1).* (samples-emu2).^2)/N2);
    end
    fprintf('mu1=%.2f mu2=%.2f s1=%.2f s2=%.2f p1=%.2f p2=%.2f\n', emu1, emu2, es1, es2,ep1, ep2);


end

