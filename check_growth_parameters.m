function check_growth_parameters()
% Check the stronger constrained-mode activation in Example 3.
% Both schedules use the same physical plant and the same stochastic paths.
folder = fileparts(mfilename('fullpath'));
p = bw_parameters();
p.A(2) = 1.25;
paths = bw_paths(p);
H = [0.25; 0.45];
q = zeros(2,1);
residuals = zeros(2,6);
generator = zeros(2,1);
clockGenerator = zeros(2,1);
terminal = zeros(2,1);
integrated = zeros(2,1);
pulseCount = zeros(2,1);
for k = 1:numel(H)
    p.H = H(k);
    p.q = exp(-p.chi*p.H)/(expm1(p.ell*p.H)-p.varpi);
    q(k) = p.q;
    cert = bw_certificate(p,false);
    residuals(k,:) = cert.residual.';
    generator(k) = cert.initialGenerator;
    clockGenerator(k) = cert.initialClockGenerator;
    run = simulate_bounded_wait(p,paths,'stic');
    assert(all(run.maxAge <= p.H+1e-10));
    terminal(k) = mean(run.terminal);
    integrated(k) = mean(run.integral);
    pulseCount(k) = mean(run.impulses);
end
% Positive physical drift at the constant initial history rules out the
% earlier common-storage zero-input dissipation estimate at both waits.
assert(all(generator > 0) && all(clockGenerator < 0));
assert(all(residuals(1,:) < 0));
assert(residuals(2,2) > 0 && all(residuals(2,[1,3,4,5,6]) < 0));
result = table(repmat(p.A(2),size(H)),H,q,generator,clockGenerator, ...
    residuals(:,1),residuals(:,2),residuals(:,3),residuals(:,4), ...
    residuals(:,5),residuals(:,6),integrated,terminal,pulseCount, ...
    'VariableNames',{'A2','WaitingBound','TriggerFloor','InitialGeneratorQ', ...
    'InitialGeneratorW','FullP','FullILower','FullIUpper','DelayP', ...
    'DelayI','Impulse','MeanIntegratedError','MeanTerminalError','MeanImpulses'});
writetable(result,fullfile(folder,'bw_growth_parameters.csv'));
disp(result);
end
