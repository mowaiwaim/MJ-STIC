function cert = verify_liu_static_certificate(p)
% Zero-input certificate for the attack-free, gated Liu static rule in Example 3.
% Q = sum_i p_i e_i^2. In mode 1 all weights are one, and in mode 2
% p_i = exp(-age_i). Each local clock resets on its node's pulse. A mode
% change resets all clocks, and the two-state generator is included below.
% Add R*integral exp(-lambda*(t-s))*||e(s)||^2 ds over the common delay.
lambda=.02; R=.016; nu=1;
G=p.Hgraph; I=eye(size(G)); N=numel(p.x0);
assert(p.A(2)==1 && p.A(1)>0 && isequal(G,G') && p.gSlope>=0);
assert(all(all(G-diag(diag(G))>=0)) && all(sum(G,2)==0));
assert(all(max(abs(p.B))*p.gSlope*ones(N,1)<=.015+eps));
P=p.M\(-p.C(1)*I+p.a*G+p.bp*p.Wgraph-p.pin);
normInv=1/min(eig(p.M));
a1=max(eig(P+P'))+2*p.A(1)*normInv+p.B(1)*p.gSlope*normInv+p.noise^2;
a2=-2*p.C(2)+2*p.A(2)+p.B(2)*p.gSlope+p.noise^2;
minWeight=exp(-nu*p.H);
% Sum over symmetric edges 2a h_ij p_i e_i(e_j-e_i). Young's
% inequality bounds the coefficient of each e_i^2 by a*sum_j h_ij(p_j-p_i).
degree=max(-diag(G));
endpoints=[minWeight,1];
flow1=a1+lambda+R;
flow2=max(endpoints*(a2+lambda-nu)+(p.Gamma(2,1)+p.a*degree)*(1-endpoints)+R);
delay=max([p.B(1)*normInv,p.B(2)])*p.gSlope-exp(-lambda*p.delay)*R;
jumpRatio=(1+p.mu)^2/minWeight;
c=norm(p.x0)^2*(1+R*(-expm1(-lambda*p.delay))/lambda)/minWeight;
assert(flow1<0 && flow2<0 && delay<0 && jumpRatio<1);
cert=table(flow1,flow2,delay,jumpRatio,lambda,c);
end
