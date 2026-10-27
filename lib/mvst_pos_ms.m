function res = mvst_pos_ms(sys, opts)
%MVST_POS_MS  Second-moment (mean-square) certificates for STOCHASTIC positive
%   systems (Chapter "Linear Positive Systems", stochastic layer): the
%   covariance-level SDP certificates that see the Brownian and Poisson noise,
%   which the mean-level copositive-LP layer (MVST_POS_STABILITY / _GAIN) is
%   blind to.  A single router over the doubly-nonnegative / semidefinite-plus-
%   nonnegative cone certificates and the diagonal quadratic certificates.
%
%   RES = MVST_POS_MS(SYS, OPTS).  OPTS.what selects the property, OPTS.cert the
%   certificate family:
%     .what : 'stability' (default) | 'tracegain' | 'peakgain'
%     .cert : 'spn' (default) | 'diag'
%
%   COVARIANCE MODEL (SYS).  Flow F(X)=He[A0 X]+sum_l Ac_l X Ac_l' with A0
%   Metzler and channels Ac_l>=0 (Brownian), plus optional Poisson jump channels
%   Jf_r>=0 at rates lamf_r contributing sum_r lamf_r (Jf_r X Jf_r' - X) to the
%   flow (prop:pos:classstoch), and optional scheduled impulsive resets
%   X(t_k^+)=sum_r J_r X(t_k^-) J_r' (J_r>=0).  The flow adjoint is
%       F*(P)=A0'P+P A0+sum_l Ac_l'P Ac_l+sum_r lamf_r (Jf_r'P Jf_r - P),
%   and the jump residual D_0(P+,P-)=sum_r J_r'P+ J_r - P- (eq:pos:dnndual and
%   the stability-chapter C_0/D_0 in the covariance metric).
%
%   'spn' CERTIFICATES (SDP + entrywise LP, EXACT cone encoding).  SPN membership
%   M in SPN_n is encoded as M = M1 + M2 with M1>=0 (PSD, SDP) and M2>=0
%   entrywise, M2 symmetric -- exactly def:pos:cones, no relaxation.  DNN
%   membership X in DNN_n is X>=0 (PSD) AND X>=0 entrywise.
%     stability  th:pos:cop: P-eps1 I in SPN, -C_0(P)-a3 I in SPN,
%                -D_0(P+,P-)-a4 I in SPN, with C_0(P)=F*(P) (constant P).
%     tracegain  prop:pos:dnntrace: the four SPN rows (eq:pos:dnntracecond),
%                minimize gamma.  Channels H_c*(I)=Cc'Cc, G_c*(P)=Ec'P Ec,
%                D_c*(I)=Fc'Fc (and the jump analogues Cd,Ed,Fd).
%     peakgain   prop:pos:dnnpeak: DNN-valued tube Xi (Xi>=0 and Xi>=0 entrywise),
%                dot/flow F(Xi)+G_c(I)<=0, output H_c(Xi)+D_c(I)<=gamma I,
%                minimize gamma.
%
%   'diag' CERTIFICATES (plain LMIs on a diagonal certificate).
%     stability  prop:pos:diag / prop:pos:dnndiag: diagonal D>0 with
%                F*(D) <= -eps I and (jump) sum_r J_r'D J_r - D <= -eps I.
%     (The structured-covariance dual action F*(D)=A0'D+D A0+Sigma D with
%     Sigma=diag(sum_l (Ac_l)_ii^2) is the diagonal case of the general F*.)
%
%   OPTS also: .eps (strict margin, default 1e-4), .solver (backend via
%   MVST_OPTS / MVST.OPTIONS, default 'sedumi'), .engine ('sostools' (default):
%   constant-degree SOSTOOLS programs (degree-zero matrix SOS variables) solved
%   through MVST_SOSSOLVE, so every backend including 'nosdp' works; the SPN
%   split M = M1 + M2 is encoded as EXPR - M2 >= 0 with the entrywise
%   nonnegative symmetric M2 | 'yalmip': legacy route, kept one release for
%   cross-checks, needs YALMIP, no 'nosdp').
%
%   RES FIELDS: .feas (certificate found), .cert ('spn'|'diag'), .what,
%     .P / .D (the certificate value), .gamma (trace/peak gain; [] for
%     stability), .margin (achieved slack), .solver, .info (YALMIP diagnostic).
%
%   See also MVST_POS_STABILITY, MVST_POS_GAIN, MVST_STABILITY, MVST_OPTS.

    if nargin < 2; opts = struct(); end
    what = lower(mvst_getf(opts,'what','stability'));
    cert = lower(mvst_getf(opts,'cert','spn'));
    engine = lower(mvst_getf(opts,'engine','sostools'));
    if ~any(strcmp(engine,{'sostools','yalmip'}))
        error('mvst:posms:engine','OPTS.engine must be ''sostools'' (default) or ''yalmip''.');
    end
    if strcmp(engine,'yalmip'); oo = mvst_opts(opts, true, 'yalmip');
    else; oo = mvst_opts(opts, true, 'sos'); end
    solver = oo.solver;
    eps0 = mvst_getf(opts,'eps',1e-4);

    A0 = mvst_getf(sys,'A0', mvst_getf(sys,'A',[]));
    if isempty(A0); error('mvst:posms:A0','SYS needs a flow drift .A0/.A.'); end
    n  = size(A0,1);
    Ac = tocell(mvst_getf(sys,'Ac', mvst_getf(sys,'Al',{})));   % Brownian channels
    Jf = tocell(mvst_getf(sys,'Jf',{}));                        % Poisson jump maps
    lamf = mvst_getf(sys,'lamf', ones(1,numel(Jf)));
    Jr = tocell(mvst_getf(sys,'J',{}));                         % scheduled resets

    res = struct('feas',false,'cert',cert,'what',what,'P',[],'D',[], ...
                 'gamma',[],'margin',[],'solver',solver,'info',[]);

    % Timer-dependent dwell-time families: the min/range dwell
    % certificates route to the SOS-timer backend (prop:pos:diagMDTRDT diagonal,
    % prop:pos:dnntrace tau/theta SPN rows).  The constant-certificate path below
    % is the arbitrary-dwell limit of these.
    dwl = lower(mvst_getf(opts,'dwell',''));
    if any(strcmp(dwl,{'min','range'}))
        res = mvst_pos_ms_dwell(sys, what, cert, opts);
        return;
    end

    if strcmp(engine,'sostools')
        switch what
            case 'stability'
                if strcmp(cert,'diag')
                    canon = mvst_getf(opts,'canonical',false);
                    if canon
                        res = ms_stab_diag(sys,A0,Ac,Jf,lamf,Jr,n,eps0,[],true,res);
                    else
                        res = ms_stab_diag_sos(A0,Ac,Jf,lamf,Jr,n,eps0,oo,res);
                    end
                else
                    res = ms_stab_spn_sos(A0,Ac,Jf,lamf,Jr,n,eps0,oo,res);
                end
            case 'tracegain'
                res = ms_trace_spn_sos(sys,A0,Ac,Jf,lamf,Jr,n,eps0,oo,res);
            case 'peakgain'
                res = ms_peak_spn_sos(sys,A0,Ac,Jf,lamf,Jr,n,eps0,oo,res);
            otherwise
                error('mvst:posms:what','OPTS.what must be stability|tracegain|peakgain.');
        end
        return;
    end

    yalmip('clear');
    yml = sdpsettings('solver',solver,'verbose',0);

    switch what
        case 'stability'
            if strcmp(cert,'diag')
                canon = mvst_getf(opts,'canonical',false);
                res = ms_stab_diag(sys,A0,Ac,Jf,lamf,Jr,n,eps0,yml,canon,res);
            else
                res = ms_stab_spn(sys,A0,Ac,Jf,lamf,Jr,n,eps0,yml,res);
            end
        case 'tracegain'
            res = ms_trace_spn(sys,A0,Ac,Jf,lamf,Jr,n,eps0,yml,res);
        case 'peakgain'
            res = ms_peak_spn(sys,A0,Ac,Jf,lamf,Jr,n,eps0,yml,res);
        otherwise
            error('mvst:posms:what','OPTS.what must be stability|tracegain|peakgain.');
    end
end

% ======================================================================= %
%  Operators                                                              %
% ======================================================================= %
function Y = Fstar(P, A0, Ac, Jf, lamf)
% F*(P) = A0'P + P A0 + sum_l Ac_l'P Ac_l + sum_r lamf_r (Jf_r'P Jf_r - P).
    Y = A0.'*P + P*A0;
    for l = 1:numel(Ac); Y = Y + Ac{l}.'*P*Ac{l}; end
    for r = 1:numel(Jf); Y = Y + lamf(r)*(Jf{r}.'*P*Jf{r} - P); end
end
function Y = Jstar(P, Jr)
% J*(P) = sum_r J_r'P J_r  (scheduled-reset jump adjoint); identity if none.
    if isempty(Jr); Y = P; return; end
    Y = zeros(size(P)); if isa(P,'sdpvar'); Y = P*0; end
    for r = 1:numel(Jr); Y = Y + Jr{r}.'*P*Jr{r}; end
end

% ======================================================================= %
%  SPN stability (th:pos:cop, constant certificate)                        %
% ======================================================================= %
function res = ms_stab_spn(sys,A0,Ac,Jf,lamf,Jr,n,eps0,yml,res) %#ok<INUSL>
    P = sdpvar(n,n,'symmetric');
    C = spn_in(P - eye(n), n);                       % P - I in SPN (P >= I, normalized scale)
    C = [C, spn_in(-Fstar(P,A0,Ac,Jf,lamf) - eps0*eye(n), n)];   % -C_0(P)-a3 I in SPN
    if ~isempty(Jr)
        C = [C, spn_in(-(Jstar(P,Jr) - P) - eps0*eye(n), n)];    % -D_0(P,P)-a4 I in SPN
    end
    d = optimize(C, [], yml);
    res.info = d;
    if d.problem == 0
        res.feas = true; res.P = value(P);
        res.margin = min(real(eig(value(-Fstar(res.P,A0,Ac,Jf,lamf)))));
    end
end

% ======================================================================= %
%  Diagonal stability (prop:pos:diag / prop:pos:dnndiag)                    %
% ======================================================================= %
function res = ms_stab_diag(sys,A0,Ac,Jf,lamf,Jr,n,eps0,yml,canon,res) %#ok<INUSL>
    if canon
        % Canonical diagonal certificate of prop:pos:dnndiag(b): the SPECIFIC
        % D = diag(e_i/d_i) with d=-A0^{-1}1, e=-A0^{-T}1 (A0 Metzler-Hurwitz,
        % Prop. prop:pos:diag).  The scaled row sums of F*(D) (via W=diag(d)) are
        % -(d_i+e_i)+d_i Sigma_i e_i, so this D certifies GUES on DNN as soon as
        % the smallness condition  d_i Sigma_i e_i < d_i + e_i  holds for every i
        % (Gershgorin, Lemma pos:gersh).  This is the exact closed-form verdict,
        % flipping precisely at the threshold; the channels must be diagonal
        % Ac_l=diag(sigma^(l)), with Sigma_i = sum_l (Ac_l)_ii^2.
        if max(real(eig(A0))) >= 0; res.feas = false; return; end
        dv = -A0\ones(n,1); ev = -A0.'\ones(n,1);
        Sig = zeros(n,1); for l = 1:numel(Ac); Sig = Sig + diag(Ac{l}).^2; end
        rowsum = dv.*Sig.*ev - (dv + ev);        % = -(d_i+e_i)+d_i Sigma_i e_i
        ok = all(rowsum < -eps0);
        res.feas = ok;
        if ok; res.D = diag(ev./dv); res.margin = -max(rowsum); end
        return;
    end
    d = sdpvar(n,1); D = diag(d);
    C = [d >= 1];                                    % D >= I (normalized, D > 0)
    C = [C, Fstar(D,A0,Ac,Jf,lamf) <= -eps0*eye(n)]; % F*(D) <= -eps I
    if ~isempty(Jr)
        C = [C, Jstar(D,Jr) - D <= -eps0*eye(n)];    % sum_r J_r'D J_r - D <= -eps I
    end
    dg = optimize(C, [], yml);
    res.info = dg;
    if dg.problem == 0
        res.feas = true; res.D = diag(value(d));
        res.margin = min(real(eig(value(-Fstar(diag(value(d)),A0,Ac,Jf,lamf)))));
    end
end

% ======================================================================= %
%  SPN trace gain (prop:pos:dnntrace, constant certificate)                %
% ======================================================================= %
function res = ms_trace_spn(sys,A0,Ac,Jf,lamf,Jr,n,eps0,yml,res) %#ok<INUSL>
    Ec = mvst_getf(sys,'Ec', mvst_getf(sys,'E',zeros(n,0)));
    Cc = mvst_getf(sys,'Cc', mvst_getf(sys,'C',zeros(0,n)));
    Fc = mvst_getf(sys,'Fc', zeros(size(Cc,1), size(Ec,2)));
    P = sdpvar(n,n,'symmetric'); gam = sdpvar(1,1);
    mc = size(Ec,2);
    C = spn_in(P - eps0*eye(n), n);                          % P in SPN (bounded below)
    C = [C, spn_in(-(Fstar(P,A0,Ac,Jf,lamf) + Cc.'*Cc), n)];        % row 1: H_c*(I)=Cc'Cc
    C = [C, spn_in(gam*eye(mc) - (Ec.'*P*Ec + Fc.'*Fc), mc)];       % row 2: G_c*(P)+D_c*(I)
    if ~isempty(Jr)
        Ed = mvst_getf(sys,'Ed', zeros(n,0)); Cd = mvst_getf(sys,'Cd', Cc);
        Fd = mvst_getf(sys,'Fd', zeros(size(Cd,1),size(Ed,2))); md = size(Ed,2);
        C = [C, spn_in(-(Jstar(P,Jr) - P + Cd.'*Cd), n)];          % row 3
        if md > 0
            C = [C, spn_in(gam*eye(md) - (Ed.'*P*Ed + Fd.'*Fd), md)]; % row 4
        end
    end
    d = optimize(C, gam, yml);
    res.info = d;
    if d.problem == 0
        res.feas = true; res.P = value(P); res.gamma = value(gam);
    end
end

% ======================================================================= %
%  DNN peak gain (prop:pos:dnnpeak, constant DNN tube)                     %
% ======================================================================= %
function res = ms_peak_spn(sys,A0,Ac,Jf,lamf,Jr,n,eps0,yml,res) %#ok<INUSL>
    Ec = mvst_getf(sys,'Ec', mvst_getf(sys,'E',zeros(n,0)));
    Cc = mvst_getf(sys,'Cc', mvst_getf(sys,'C',zeros(0,n)));
    Fc = mvst_getf(sys,'Fc', zeros(size(Cc,1), size(Ec,2)));
    Xi = sdpvar(n,n,'symmetric'); gam = sdpvar(1,1);
    qc = size(Cc,1); mc = size(Ec,2);
    % DNN tube: Xi in DNN (PSD and entrywise nonnegative)
    C = [Xi >= eps0*eye(n), Xi(:) >= 0];
    % flow: F(Xi) + G_c(I) <= 0  (primal covariance tube grows no faster)
    FXi = A0*Xi + Xi*A0.';
    for l = 1:numel(Ac); FXi = FXi + Ac{l}*Xi*Ac{l}.'; end
    for r = 1:numel(Jf); FXi = FXi + lamf(r)*(Jf{r}*Xi*Jf{r}.' - Xi); end
    C = [C, FXi + Ec*Ec.' <= 0];                              % G_c(I)=Ec Ec'
    C = [C, Cc*Xi*Cc.' + Fc*Fc.' <= gam*eye(qc)];             % H_c(Xi)+D_c(I) <= gamma I
    if ~isempty(Jr)
        JXi = zeros(n); for r=1:numel(Jr); JXi = JXi + Jr{r}*Xi*Jr{r}.'; end
        C = [C, JXi - Xi <= 0];                               % jump non-expansive on the tube
    end
    d = optimize(C, gam, yml);
    res.info = d;
    if d.problem == 0
        res.feas = true; res.P = value(Xi); res.gamma = value(gam);
    end
end

% ======================================================================= %
%  SOSTOOLS engine (degree-zero matrix SOS variables, MVST_SOSSOLVE)       %
% ======================================================================= %
function [prog,ok] = sos_run(prog, solver, obj)
    prog = psd_flush(prog);
    if nargin >= 3; prog = sossetobj(prog, obj); end
    o = mvst_sosopts(solver);
    [prog,info] = mvst_sossolve(prog, o);
    ok = (info.pinf == 0 && info.dinf == 0 && info.numerr < 2);
end
function [prog,M] = sos_sym(prog, th, n)
    [prog,M] = sospolymatrixvar(prog, monomials(th,0), [n n], 'symmetric');
end
function prog = spn_in_sos(prog, th, EXPR, n)
% EXPR in SPN_n  <=>  EXPR - M2 >= 0 (PSD) for some symmetric M2 >= 0 entrywise.
    [prog,M2] = sos_sym(prog, th, n);
    for i = 1:n
        for j = i:n; prog = psd_ineq(prog, M2(i,j)); end    % (M2)_ij >= 0
    end
    prog = psd_ineq(prog, EXPR - M2);
end

function res = ms_stab_spn_sos(A0,Ac,Jf,lamf,Jr,n,eps0,oo,res)
    pvar th; prog = sosprogram(th);
    [prog,P] = sos_sym(prog, th, n);
    prog = spn_in_sos(prog, th, P - eye(n), n);
    prog = spn_in_sos(prog, th, -Fstar(P,A0,Ac,Jf,lamf) - eps0*eye(n), n);
    if ~isempty(Jr)
        prog = spn_in_sos(prog, th, -(Jstar(P,Jr) - P) - eps0*eye(n), n);
    end
    [prog,ok] = sos_run(prog, oo);
    res.info = struct('engine','sostools','ok',ok);
    if ok
        res.feas = true; res.P = double(sosgetsol(prog,P));
        res.margin = min(real(eig(-Fstar(res.P,A0,Ac,Jf,lamf))));
    end
end

function res = ms_stab_diag_sos(A0,Ac,Jf,lamf,Jr,n,eps0,oo,res)
    pvar th; prog = sosprogram(th);
    [prog,D] = sos_sym(prog, th, n);
    for i = 1:n
        prog = psd_ineq(prog, D(i,i) - 1);                  % d_i >= 1
        for j = i+1:n; prog = soseq(prog, D(i,j)); end           % D diagonal
    end
    prog = psd_ineq(prog, -Fstar(D,A0,Ac,Jf,lamf) - eps0*eye(n));   % F*(D) <= -eps I
    if ~isempty(Jr)
        prog = psd_ineq(prog, -(Jstar(D,Jr) - D) - eps0*eye(n));    % J*(D)-D <= -eps I
    end
    [prog,ok] = sos_run(prog, oo);
    res.info = struct('engine','sostools','ok',ok);
    if ok
        Dv = double(sosgetsol(prog,D)); Dv = diag(diag(Dv));
        res.feas = true; res.D = Dv;
        res.margin = min(real(eig(-Fstar(Dv,A0,Ac,Jf,lamf))));
    end
end

function res = ms_trace_spn_sos(sys,A0,Ac,Jf,lamf,Jr,n,eps0,oo,res)
    Ec = mvst_getf(sys,'Ec', mvst_getf(sys,'E',zeros(n,0)));
    Cc = mvst_getf(sys,'Cc', mvst_getf(sys,'C',zeros(0,n)));
    Fc = mvst_getf(sys,'Fc', zeros(size(Cc,1), size(Ec,2)));
    mc = size(Ec,2);
    pvar th; prog = sosprogram(th); gam = dpvar('gam'); prog = sosdecvar(prog,gam);
    [prog,P] = sos_sym(prog, th, n);
    prog = spn_in_sos(prog, th, P - eps0*eye(n), n);
    prog = spn_in_sos(prog, th, -(Fstar(P,A0,Ac,Jf,lamf) + Cc.'*Cc), n);
    prog = spn_in_sos(prog, th, gam*eye(mc) - (Ec.'*P*Ec + Fc.'*Fc), mc);
    if ~isempty(Jr)
        Ed = mvst_getf(sys,'Ed', zeros(n,0)); Cd = mvst_getf(sys,'Cd', Cc);
        Fd = mvst_getf(sys,'Fd', zeros(size(Cd,1),size(Ed,2))); md = size(Ed,2);
        prog = spn_in_sos(prog, th, -(Jstar(P,Jr) - P + Cd.'*Cd), n);
        if md > 0
            prog = spn_in_sos(prog, th, gam*eye(md) - (Ed.'*P*Ed + Fd.'*Fd), md);
        end
    end
    [prog,ok] = sos_run(prog, oo, gam);
    res.info = struct('engine','sostools','ok',ok);
    if ok
        res.feas = true; res.P = double(sosgetsol(prog,P)); res.gamma = double(sosgetsol(prog,gam));
    end
end

function res = ms_peak_spn_sos(sys,A0,Ac,Jf,lamf,Jr,n,eps0,oo,res)
    Ec = mvst_getf(sys,'Ec', mvst_getf(sys,'E',zeros(n,0)));
    Cc = mvst_getf(sys,'Cc', mvst_getf(sys,'C',zeros(0,n)));
    Fc = mvst_getf(sys,'Fc', zeros(size(Cc,1), size(Ec,2)));
    qc = size(Cc,1);
    pvar th; prog = sosprogram(th); gam = dpvar('gam'); prog = sosdecvar(prog,gam);
    [prog,Xi] = sos_sym(prog, th, n);
    prog = psd_ineq(prog, Xi - eps0*eye(n));                 % Xi PSD (>= eps I)
    for i = 1:n
        for j = i:n; prog = psd_ineq(prog, Xi(i,j)); end     % entrywise Xi >= 0
    end
    FXi = A0*Xi + Xi*A0.';
    for l = 1:numel(Ac); FXi = FXi + Ac{l}*Xi*Ac{l}.'; end
    for r = 1:numel(Jf); FXi = FXi + lamf(r)*(Jf{r}*Xi*Jf{r}.' - Xi); end
    prog = psd_ineq(prog, -(FXi + Ec*Ec.'));                 % F(Xi)+G_c(I) <= 0
    prog = psd_ineq(prog, gam*eye(qc) - (Cc*Xi*Cc.' + Fc*Fc.'));
    if ~isempty(Jr)
        JXi = zeros(n); for r=1:numel(Jr); JXi = JXi + Jr{r}*Xi*Jr{r}.'; end
        prog = psd_ineq(prog, -(JXi - Xi));                  % jump non-expansive on the tube
    end
    [prog,ok] = sos_run(prog, oo, gam);
    res.info = struct('engine','sostools','ok',ok);
    if ok
        res.feas = true; res.P = double(sosgetsol(prog,Xi)); res.gamma = double(sosgetsol(prog,gam));
    end
end

function prog = psd_ineq(prog, EXPR)
% EXPR >= 0 as a matrix SOS constraint.  1x1 constraints are collected in
% prog.mvst_pending and flushed as ONE diagonal block (see PSD_FLUSH): programs
% programs whose Gram blocks are ALL 1x1 crash SeDuMi's pretransfo (still the case;
% the former SOSNOSDP crash on mixed 1x1/larger blocks, nsd_gram>apply_A, is fixed in
% the vendored sosnosdp, so this merge is now needed for SeDuMi only).
    if size(EXPR,1) == 1
        if ~isfield(prog,'mvst_pending'); prog.mvst_pending = {}; end
        prog.mvst_pending{end+1} = EXPR;
    else
        prog = sosmatrixineq(prog, EXPR);
    end
end

function prog = psd_flush(prog)
    if ~isfield(prog,'mvst_pending') || isempty(prog.mvst_pending); return; end
    L = prog.mvst_pending; k = numel(L); if k == 1; k = 2; end
    D = 0 * L{1} * zeros(k);
    for i = 1:numel(L); Ei = zeros(k); Ei(i,i) = 1; D = D + L{i}*Ei; end
    if numel(L) == 1; D = D + 1*[0 0;0 1]; end
    prog = rmfield(prog,'mvst_pending');
    prog = sosmatrixineq(prog, D);
end

% ======================================================================= %
%  Cone-membership encodings (EXACT, per def:pos:cones)                    %
% ======================================================================= %
function C = spn_in(EXPR, n)
% EXPR in SPN_n  <=>  EXPR = M1 + M2, M1 >= 0 (PSD), M2 symmetric entrywise >= 0.
    M1 = sdpvar(n,n,'symmetric'); M2 = sdpvar(n,n,'symmetric');
    C = [M1 >= 0, M2(:) >= 0, EXPR == M1 + M2];
end

% ======================================================================= %
function c = tocell(X)
    if iscell(X); c = X; elseif isempty(X); c = {};
    else; p = size(X,3); c = cell(1,p); for k=1:p; c{k}=X(:,:,k); end; end
end
