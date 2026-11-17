function paths = bw_paths(p)
% Exact exponential CTMC residence times, separate Brownian and CTMC streams.
n = round(p.T/p.dt);
assert(abs(n*p.dt-p.T)<1e-10);
rsNoise = RandStream('mt19937ar','Seed',p.seed);
rsMode = RandStream('mt19937ar','Seed',p.seed+1);
paths.dB = sqrt(p.dt)*randn(rsNoise,n,p.samples);
paths.switches = cell(p.samples,1);
for trial = 1:p.samples
    times = [];
    t = 0;
    mode = 2;
    while true
        t = t - log(rand(rsMode))/(-p.Gamma(mode,mode));
        if t >= p.T, break; end
        times(end+1) = t; %#ok<AGROW>
        mode = 3-mode;
    end
    paths.switches{trial}=times;
end
end
