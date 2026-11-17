function cert = bw_certificate(p,requireFeasible)
% Check the exponential output-strict passivity conditions for arbitrary input.
% The delay block increases with the weight; check it at c_r in each mode.
% Supply false to inspect the residual of a trial outside this certificate.
if nargin < 2, requireFeasible = true; end
N=numel(p.x0);
I=eye(N);
Pi=I-ones(N)/N;
M=p.M;
rho=p.rho;
Qfp=2*p.A(1);
Qfi=2*p.A(2)*norm(M,2);
theta=p.gSlope^2;
minWeight=p.c(2)*exp(-p.nu*p.H);
currentP=(-2*p.C(1)+Qfp+p.B(1)^2+rho*p.R1)*I ...
    +2*(p.a*p.Hgraph+p.bp*p.Wgraph)-2*p.pin;
linearI=M*(-p.C(2)*I+p.a*p.Hgraph);
currentI=linearI+linearI'+Qfi*I+0.5*p.B(2)^2*M*M+rho*p.R1*I;

% Each row specifies a mode and an endpoint of its timer-weight interval.
modeAtEndpoint=[1;2;2];
weights=[p.c(1);minWeight;p.c(2)];
current={currentP,currentI};
input={Pi,M*Pi};
nu=[0 p.nu];
E=p.E1*I;
D=p.E2*I;
V=-(D'*p.Z+p.Z'*D)+2*D'*p.Lambda2*D;
worst=@(A) max(eig((A+A')/2));
fullBlock=zeros(3,1);
for k=1:numel(weights)
    r=modeAtEndpoint(k);
    w=weights(k);
    j=3-r;
    X=w*(current{r}+(p.lambda-nu(r))*M)+p.R ...
        +p.Gamma(r,j)*(p.c(j)-w)*M+2*E'*p.Lambda2*E;
    U=w*input{r}-E'*p.Z+2*E'*p.Lambda2*D;
    fullBlock(k)=worst([X,U;U',V]);
end
delayP=p.c(1)*(rho*p.R2+theta)*I-exp(-p.lambda*p.delay)*p.R;
delayI=p.c(2)*(rho*p.R2+2*theta)*I-exp(-p.lambda*p.delay)*p.R;
jump=(1+p.mu)^2*exp(p.nu*p.H)-1;
cert.residual=[fullBlock;worst(delayP);worst(delayI);jump];
cert.Qfp=Qfp;
cert.Qfi=Qfi;
cert.minM=min(eig(M));
cert.maxM=worst(M);

% For the subsequent zero-input experiment the leader remains at zero.
Q0=p.x0'*M*p.x0;
drift=(-p.C(2)*I+p.a*p.Hgraph)*p.x0+p.A(2)*tanh(p.x0) ...
    +p.B(2)*p.gSlope*tanh(p.x0);
cert.initialGenerator=2*p.x0'*M*drift+p.noise^2*Q0;
% At the initial constant history the delay term's derivative is zero.
% Mode 2 starts at timer age zero; the two reset weights are both one.
cert.initialClockGenerator=p.c(2)*(cert.initialGenerator-p.nu*Q0) ...
    +p.Gamma(2,1)*(p.c(1)-p.c(2))*Q0;
initialDelay=p.x0'*p.R*p.x0*(-expm1(-p.lambda*p.delay))/p.lambda;
cert.boundConstant=(p.c(2)*Q0+initialDelay)/(min(p.c(1),minWeight)*cert.minM);
if requireFeasible
    assert(all(cert.residual<0) && cert.initialGenerator>0);
end
assert(cert.minM>0 && cert.maxM<=rho && min(eig(p.Lambda2))>0);
assert(p.q>0 && log1p(p.varpi)/p.ell<p.H && min(eig(p.R))>=0);
assert(abs(p.q*(expm1(p.ell*p.H)-p.varpi)-exp(-p.chi*p.H))<1e-13);
end
