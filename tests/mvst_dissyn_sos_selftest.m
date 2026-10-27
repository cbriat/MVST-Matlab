function res = mvst_dissyn_sos_selftest()
%MVST_DISSYN_SOS_SELFTEST  Validation of MVST_DISSYN_SOS by closed-loop
%   re-analysis.  Each case SYNTHESIZES a state-feedback gain with
%   MVST_DISSYN_SOS, closes the loop with the recovered gain, and RE-ANALYSES
%   the closed loop with the exact backend MVST_L2GAIN_CT and/or the analysis
%   routine MVST_DISS_SOS.  The achieved dissipativity level must not exceed
%   the certified synthesis level (and, for the CT LTI L2-gain relaxation which
%   is exact, must agree with it to solver tolerance).  Returns a struct array
%   (.name,.expected,.got,.ok) and prints a table.
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
addpath(root); addpath(fullfile(root,'lib'));
if exist('sossolve','file') ~= 2 || exist('sedumi','file') == 0
    fprintf('%s: SKIPPED (SOSTOOLS + SeDuMi required; run mvst_setup).\n', mfilename);
    res = [];
    return
end

R = struct('name',{},'expected',{},'got',{},'ok',{});

% ---- (1) CT arbitrary L2-gain synthesis, exact cross-check --------------
% Plant  xdot = A0 x + B0 u + Bw w,  z = [x1 ; u].  The constant-W synthesis
% relaxation is EXACT for the CT LTI H-infinity state-feedback problem, so the
% H-infinity norm of the synthesized closed loop equals the certified level.
% Two disturbance channels (m = 2) match the two outputs (pz = 2), so the exact
% Hamiltonian backend MVST_L2GAIN_CT (which assumes square output/disturbance)
% applies directly.
A0=[0 1;-2 -1]; B0=[0;1]; Bw=eye(2);
C=[1 0;0 0]; Du=[0;1]; Dw=zeros(2,2);
sys1 = struct('A0',A0,'B0',B0,'Bw',Bw,'C',C,'Du',Du,'Dw',Dw);
r1 = mvst_dissyn_sos(sys1, struct('type','l2gain'), struct('type','arbitrary'), 'gamma');
K1  = r1.K(0);
Acl = A0 + B0*K1;   Ccl = C + Du*K1;
gExact = mvst_l2gain_ct(Acl, Bw, Ccl, Dw);
R(end+1) = rowub('1 CT L2-gain synth (exact)', gExact, r1.gammastar, 4e-2);

% ---- (2) recovered central gain stabilizes the closed loop --------------
R(end+1) = rowb('2 central gain stabilizes', true, max(real(eig(Acl))) < 0);

% ---- (3) CT arbitrary passivity synthesis, closed-loop re-analysis ------
% Oscillator with collocated (velocity) output; synthesize a passifying gain
% then confirm the analysis routine certifies the closed loop passive.
A0=[0 1;-1 0]; B0=[0;1]; Bw=[0;1]; C=[0 1]; Du=0; Dw=0.5;
sys3 = struct('A0',A0,'B0',B0,'Bw',Bw,'C',C,'Du',Du,'Dw',Dw);
r3 = mvst_dissyn_sos(sys3, struct('type','passivity'), struct('type','arbitrary'), 'feas');
K3  = r3.K(0);
cl3 = struct('A0',A0+B0*K3,'Bw',Bw,'C',C+Du*K3,'Dw',Dw);
a3  = mvst_diss_sos(cl3, struct('type','passivity'), struct('type','arbitrary'), 'feas');
R(end+1) = rowb('3 passivity synth (feas+reanal)', true, r3.feas && a3.feas);

% ---- (4) periodic sampled-output L2-gain synthesis (constant gain) ------
% Flow carries no output (pure disturbance credit); the output is sampled at
% the periodic jump (CJ).  Constant certificate (deg 0) so the closed loop has
% constant matrices and can be re-analysed exactly by the periodic analysis
% DLMIs.  The analysed achieved gain must not exceed the certified level.
A0=[0 1;-1 -1]; B0=[0;1]; Bw=[0;1]; CJ=[1 0]; Tbar=0.3;
sys4 = struct('A0',A0,'B0',B0,'Bw',Bw,'J',eye(2),'CJ',CJ);
r4 = mvst_dissyn_sos(sys4, struct('type','l2gain'), ...
                     struct('type','periodic','Tbar',Tbar,'deg',0), 'gamma');
K4  = r4.K(0);
cl4 = struct('A0',A0+B0*K4,'Bw',Bw,'J',eye(2),'CJ',CJ);
g4a = mvst_diss_sos(cl4, struct('type','l2gain'), ...
                    struct('type','periodic','Tbar',Tbar,'deg',4), 'gamma');
R(end+1) = rowle('4 periodic synth (achieved<=cert)', r4.gammastar, g4a.gamma, 4e-2);

% ---- (5) timer-polynomial refinement (validates the interval SOS path) --
% A degree-4 timer certificate has strictly more freedom than a constant one,
% so its certified level cannot exceed the deg-0 level: exercises the matrix-SOS
% multiplier machinery on [0,Tbar] and the refinement monotonicity.
r4d4 = mvst_dissyn_sos(sys4, struct('type','l2gain'), ...
                       struct('type','periodic','Tbar',Tbar,'deg',4), 'gamma');
R(end+1) = rowle('5 periodic deg4 <= deg0', r4.gammastar, r4d4.gammastar, 4e-2);

% ---- (6) minimum dwell-time L2-gain synthesis, closed-loop re-analysis --
% Flow output C x with feedthrough Dw (gamma floor = norm(Dw)) and a genuine
% contraction reset J at the dwell.  Constant certificate (deg 0) so the closed
% loop is constant and re-analysable by the minimum-dwell analysis DLMIs.
A0=[0 1;-1 -1]; B0=[0;1]; Bw=[0;1]; C=[1 0]; Dw=0.2; J=0.9*eye(2); Tmin=0.4;
sys6 = struct('A0',A0,'B0',B0,'Bw',Bw,'C',C,'Dw',Dw,'J',J);
mdt  = struct('type','minimum','Tbar',Tmin,'deg',0);
r6 = mvst_dissyn_sos(sys6, struct('type','l2gain'), mdt, 'gamma');
K6  = r6.K(0);
cl6 = struct('A0',A0+B0*K6,'Bw',Bw,'C',C,'Dw',Dw,'J',J);
g6a = mvst_diss_sos(cl6, struct('type','l2gain'), ...
                    struct('type','minimum','Tbar',Tmin,'deg',4), 'gamma');
R(end+1) = rowle('6 min-dwell synth (achieved<=cert)', r6.gammastar, g6a.gamma, 5e-2);

% ---- (7-11) five-seed exactness battery -------------------------
% Randomized CT arbitrary L2-gain state-feedback instances (n=3, two outputs,
% two disturbances).  The constant-W relaxation is EXACT for CT LTI H-infinity
% state feedback, so on the well-conditioned seeds 1,2,3,5 the certified level
% must match the closed-loop H-infinity norm (mvst_l2gain_ct) to 1e-3.  Seed 4
% is the adversarial instance: before the post-solve validation gate
% mvst_dissyn_sos returned feas=true with an INDEFINITE certificate, an
% unsatisfied flow LMI, and a DESTABILIZING gain (abscissa +4873).  The gate
% must make that impossible: the routine must either return feas=false at the
% previously-lying level or a genuinely VALID result (W>0, stable closed loop,
% achieved <= certified*1.05).  Every seed row also re-checks validity.
for s = [1 2 3 4 5]
    rng(s);
    A0=randn(3)-2*eye(3); B0=randn(3,1); Bw=randn(3,2);
    C=randn(2,3); Du=[0;1]; Dw=zeros(2,2);
    syss = struct('A0',A0,'B0',B0,'Bw',Bw,'C',C,'Du',Du,'Dw',Dw);
    rs = mvst_dissyn_sos(syss, struct('type','l2gain'), ...
                         struct('type','arbitrary'), 'gamma');
    if rs.feas
        Ks = rs.K(0); Acl = A0+B0*Ks; Ccl = C+Du*Ks; Ws = rs.W(0);
        pdOK   = min(eig((Ws+Ws')/2)) > 0;          % certificate positive def.
        stabOK = max(real(eig(Acl))) < 0;           % closed loop stable
        ach    = mvst_l2gain_ct(Acl, Bw, Ccl, Dw);  % exact re-analysis
        boundOK = ach <= rs.gammastar*1.05;         % achieved within 5% of cert
        valid  = pdOK && stabOK && boundOK;
        cert   = rs.gammastar;
    else
        valid = true; ach = NaN; cert = NaN;        % legitimate infeasible is OK
    end
    if s == 4
        % Regression: NO LIE.  Either feasible-and-valid, or infeasible.
        R(end+1) = rowseed('10 seed-4 regression (no lie)', cert, ach, valid); %#ok<AGROW>
    else
        % Exactness: valid AND certified level matches achieved to 1e-3.
        exact = rs.feas && (abs(ach-cert) <= 1e-3*max(1,cert));
        R(end+1) = rowseed(sprintf('%d seed-%d exact (|d|<=1e-3)', 6+s, s), ...
                           cert, ach, valid && exact); %#ok<AGROW>
    end
end

% ---- print --------------------------------------------------------------
fprintf('\n==== mvst_dissyn_sos self-test (synthesis vs closed-loop re-analysis) ====\n');
fprintf('%-34s %12s %12s %6s\n','case','reference','synth/got','ok');
for i=1:numel(R)
    if isnan(R(i).expected)
        fprintf('%-34s %12s %12.6f %6d\n', R(i).name, '(bool)', R(i).got, R(i).ok);
    else
        fprintf('%-34s %12.6f %12.6f %6d\n', R(i).name, R(i).expected, R(i).got, R(i).ok);
    end
end
fprintf('ALL PASS: %d\n', all([R.ok]));
res = R;
end

function s = rowub(name, gref, ggot, tol)
% exact-relaxation row: achieved gain gref must be <= certified ggot (up to
% tol) AND agree with it to tol (the CT LTI relaxation is exact).
s.name=name; s.expected=gref; s.got=ggot;
s.ok = (gref <= ggot + tol*max(1,abs(ggot))) && ...
       (abs(ggot-gref) <= tol*max(1,abs(gref)));
end
function s = rowle(name, gcert, gachieved, tol)
% upper-bound row: analysed achieved level must not exceed the certified level.
s.name=name; s.expected=gcert; s.got=gachieved;
s.ok = gachieved <= gcert*(1+tol) + tol;
end
function s = rowb(name, expectedBool, gotBool)
s.name=name; s.expected=NaN; s.got=double(gotBool);
s.ok = (logical(gotBool)==logical(expectedBool));
end
function s = rowseed(name, cert, achieved, okFlag)
% seed-battery row: 'reference' column shows the certified level, 'got' the
% closed-loop achieved gain; ok is decided by the caller (validity/exactness).
s.name=name; s.expected=cert; s.got=achieved; s.ok = logical(okFlag);
end
