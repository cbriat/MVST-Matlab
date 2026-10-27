function res = mvst_pos_ms_dwell(sys, what, cert, opts)
%MVST_POS_MS_DWELL  Timer-dependent (minimum / range dwell-time) second-moment
%   certificates for stochastic positive systems, the SOS-timer
%   backend of MVST_POS_MS.  Extends the constant-certificate SPN/diagonal
%   certificates to the min/range dwell families of prop:pos:diagMDTRDT (diagonal
%   quadratic), prop:pos:dnntrace (SPN trace gain, tau/theta rows), and
%   prop:pos:dnnpeak + th:pos:peakgainthm (DNN peak gain, eq:pos:peakcond tube).
%
%   RES = MVST_POS_MS_DWELL(SYS, WHAT, CERT, OPTS), WHAT in
%   {stability,tracegain,peakgain}, CERT in {diag,spn}, OPTS.dwell in {min,range}.
%
%   diag / stability (prop:pos:diagMDTRDT): diagonal D(tau)>0 with
%     dot D + F*(D) <= -eps I   on (0,T),   J' D(0+) J - D(theta) <= -eps I
%   (theta=Tmin for min, all theta in [Tmin,Tmax] for range), and for minimum
%   dwell-time the frozen tail  F*(D(Tmin)) <= -eps I.
%   F*(D)=A0'D+D A0+sum_l Ac_l'D Ac_l+sum_r lamf_r(Jf_r'D Jf_r - D).
%
%   spn / tracegain (prop:pos:dnntrace): P(tau) in SPN with, for all tau,theta,
%     -(dot P + F*(P) + Cc'Cc) in SPN,        gamma I - (Ec'P Ec + Fc'Fc) in SPN,
%     -(J*(P(0+)) - P(theta) + Cd'Cd) in SPN, gamma I - (Ed'P(0+)Ed + Fd'Fd) in SPN,
%   with the frozen tail -(F*(P(Tmin))+Cc'Cc) in SPN for minimum dwell-time;
%   minimize gamma.  SPN membership of a timer matrix M(tau) is the EXACT cone
%   split  M = M1 + M2,  M1(tau) matrix-SOS (PSD on [0,T]),  M2(tau) symmetric
%   with every entry a nonnegative polynomial on [0,T] (def:pos:cones).
%
%   OPTS: .dwell, .Tmin, .Tmax, .deg (timer degree, default 4), .eps, .solver.
%   RES: .feas, .gamma (tracegain), .Pfun/.Dfun (@(tau)), .cert, .what, .dwell.
%
%   See also MVST_POS_MS, MVST_STAB_SOS, MVST_SOSPOSMATRVAR.

    oo = mvst_opts(opts, true, 'sos'); solver = oo.solver;
    A0 = mvst_getf(sys,'A0', mvst_getf(sys,'A',[]));
    n  = size(A0,1);
    Ac = tocell(mvst_getf(sys,'Ac', mvst_getf(sys,'Al',{})));
    Jf = tocell(mvst_getf(sys,'Jf',{})); lamf = mvst_getf(sys,'lamf', ones(1,numel(Jf)));
    Jr = tocell(mvst_getf(sys,'J',{}));
    deg  = mvst_getf(opts,'deg',4);  eps0 = mvst_getf(opts,'eps',1e-4);
    dwl  = lower(mvst_getf(opts,'dwell',''));
    Tmin = mvst_getf(opts,'Tmin',[]); Tmax = mvst_getf(opts,'Tmax',Tmin);
    if strcmp(dwl,'range'); T=Tmax; else; T=Tmin; end
    res = struct('feas',false,'gamma',[],'Pfun',[],'Dfun',[],'cert',cert, ...
                 'what',what,'dwell',dwl,'info',[]);

    pvar th; gt=th*(1-th);
    if strcmp(dwl,'range'); gm=(th-Tmin/T)*(1-th); else; gm=[]; end

    if strcmp(cert,'diag') && strcmp(what,'stability')
        sl=dpvar('sl'); prog=sosprogram(th); prog=sosdecvar(prog,sl);
        prog=sosineq(prog,1-sl);                                      % slack <= 1: scaling a feasible (D,sl) by c>1 stays feasible (monotonicity), so the maximized slack is unbounded without this cap
        [prog,dvec]=sospolymatrixvar(prog,monomials(th,0:deg),[n 1]);
        D=mkdiag(dvec,n); Dd=mkdiag(diff(dvec,th),n)/T;
        Fs=Fstar_ms(D,A0,Ac,Jf,lamf);
        prog=mat_ge(prog,-(Dd+Fs)-sl*eye(n),gt,th,deg,n);          % dotD+F*(D)<=-sl I
        prog=mat_ge(prog,D-eps0*eye(n),gt,th,deg,n);               % D>=eps I
        D0=subs(D,th,0); DT=subs(D,th,1);
        if ~isempty(Jr)
            prog=jump_diag(prog,D0,D,DT,Jr,sl,dwl,gm,th,deg,n);    % reset jump rows
        end
        if strcmp(dwl,'min')
            prog=sosmatrixineq(prog,-Fstar_ms(DT,A0,Ac,Jf,lamf)-sl*eye(n)); % frozen tail
        end
        prog=sossetobj(prog,-sl); o=mvst_sosopts(oo); prog=mvst_sossolve(prog, o);
        s=double(sosgetsol(prog,sl)); Dv=sosgetsol(prog,D);
        res.feas = s>1e-7; res.Dfun=@(tau) double(subs(Dv,th,min(max(tau,0),T)/T));
        res.info=struct('slack',s);

    elseif strcmp(cert,'spn') && strcmp(what,'tracegain')
        Ec=mvst_getf(sys,'Ec',mvst_getf(sys,'E',zeros(n,0)));
        Cc=mvst_getf(sys,'Cc',mvst_getf(sys,'C',zeros(0,n)));
        Fc=mvst_getf(sys,'Fc',zeros(size(Cc,1),size(Ec,2))); mc=size(Ec,2);
        g=dpvar('g'); prog=sosprogram(th); prog=sosdecvar(prog,g);
        [prog,P]=sospolymatrixvar(prog,monomials(th,0:deg),[n n],'symmetric');
        Pd=diff(P,th)/T;
        prog=spn_ge(prog,P-eps0*eye(n),gt,th,deg,n);               % P in SPN (bounded below)
        prog=spn_ge(prog,-(Pd+Fstar_ms(P,A0,Ac,Jf,lamf)+Cc.'*Cc),gt,th,deg,n);  % row1
        prog=spn_ge(prog,g*eye(mc)-(Ec.'*P*Ec+Fc.'*Fc),gt,th,deg,mc);           % row2
        P0=subs(P,th,0); PT=subs(P,th,1);
        if ~isempty(Jr)
            prog=jump_spn(prog,P0,P,PT,Jr,sys,g,dwl,gm,th,deg,n);  % rows 3,4
        end
        if strcmp(dwl,'min')
            prog=spn_pt(prog,-(Fstar_ms(PT,A0,Ac,Jf,lamf)+Cc.'*Cc),th,n);  % frozen tail
        end
        prog=sossetobj(prog,g); o=mvst_sosopts(oo); prog=mvst_sossolve(prog, o);
        gv=double(sosgetsol(prog,g)); Pv=sosgetsol(prog,P);
        res.feas=isfinite(gv); res.gamma=gv; res.Pfun=@(tau) double(subs(Pv,th,min(max(tau,0),T)/T));

    elseif strcmp(cert,'spn') && strcmp(what,'peakgain')
        % DNN peak gain (prop:pos:dnnpeak + th:pos:peakgainthm, eq:pos:peakcond):
        % DNN-valued tube Xi(tau) (PSD AND entrywise >=0) with the primal flow,
        % jump, and output rows; frozen tail eq:pos:pktail for minimum dwell-time.
        Ec=mvst_getf(sys,'Ec',mvst_getf(sys,'E',zeros(n,0)));
        Cc=mvst_getf(sys,'Cc',mvst_getf(sys,'C',zeros(0,n)));
        Fc=mvst_getf(sys,'Fc',zeros(size(Cc,1),size(Ec,2))); qc=size(Cc,1);
        Ed=mvst_getf(sys,'Ed',zeros(n,0)); Cd=mvst_getf(sys,'Cd',zeros(0,n));
        Fd=mvst_getf(sys,'Fd',zeros(size(Cd,1),size(Ed,2))); qd=size(Cd,1);
        g=dpvar('g'); prog=sosprogram(th); prog=sosdecvar(prog,g);
        [prog,Xi]=sospolymatrixvar(prog,monomials(th,0:deg),[n n],'symmetric');
        Xid=diff(Xi,th)/T;
        prog=mat_ge(prog,Xi-eps0*eye(n),gt,th,deg,n);              % Xi PSD (>= eps I)
        prog=entrywise_ge(prog,Xi,gt,th,deg,n);                   % Xi entrywise >= 0 (DNN)
        prog=mat_ge(prog,Xid-Fprimal_ms(Xi,A0,Ac,Jf,lamf)-Ec*Ec.',gt,th,deg,n); % pkflow
        prog=mat_ge(prog,g*eye(qc)-Cc*Xi*Cc.'-Fc*Fc.',gt,th,deg,qc);            % pkoutc
        Xi0=subs(Xi,th,0); XiT=subs(Xi,th,1);
        if strcmp(dwl,'range')
            prog=mat_ge(prog,Xi0-Jprimal_ms(Xi,Jr)-Ed*Ed.',gm,th,deg,n);        % pkjump, all theta
            if qd>0; prog=mat_ge(prog,g*eye(qd)-Cd*Xi*Cd.'-Fd*Fd.',gm,th,deg,qd); end % pkoutd
        else
            prog=sosmatrixineq(prog,Xi0-Jprimal_ms(XiT,Jr)-Ed*Ed.');            % pkjump, theta=Tmin
            if qd>0; prog=sosmatrixineq(prog,g*eye(qd)-Cd*XiT*Cd.'-Fd*Fd.'); end % pkoutd
            prog=sosmatrixineq(prog,-(Fprimal_ms(XiT,A0,Ac,Jf,lamf)+Ec*Ec.'));  % pktail
        end
        prog=sossetobj(prog,g); o=mvst_sosopts(oo); prog=mvst_sossolve(prog, o);
        gv=double(sosgetsol(prog,g)); Xv=sosgetsol(prog,Xi);
        res.feas=isfinite(gv); res.gamma=gv; res.Pfun=@(tau) double(subs(Xv,th,min(max(tau,0),T)/T));
    else
        error('mvst:posms:dwell','unsupported (what,cert)=(%s,%s) for dwell.',what,cert);
    end
end

% ---- operators / jump rows -------------------------------------------- %
function Y = Fstar_ms(P,A0,Ac,Jf,lamf)
    Y = A0.'*P+P*A0;
    for l=1:numel(Ac); Y=Y+Ac{l}.'*P*Ac{l}; end
    for r=1:numel(Jf); Y=Y+lamf(r)*(Jf{r}.'*P*Jf{r}-P); end
end
function Y = Jstar_ms(P,Jr)
    if isempty(Jr); Y=P; return; end
    Y=P*0; for r=1:numel(Jr); Y=Y+Jr{r}.'*P*Jr{r}; end
end
function Y = Fprimal_ms(X,A0,Ac,Jf,lamf)
% Primal covariance flow F(X)=A0 X+X A0'+sum_l Ac_l X Ac_l'+sum_r lamf_r(Jf_r X Jf_r'-X).
    Y = A0*X+X*A0.';
    for l=1:numel(Ac); Y=Y+Ac{l}*X*Ac{l}.'; end
    for r=1:numel(Jf); Y=Y+lamf(r)*(Jf{r}*X*Jf{r}.'-X); end
end
function Y = Jprimal_ms(X,Jr)
% Primal reset J(X)=sum_r J_r X J_r'.
    if isempty(Jr); Y=X*0; return; end
    Y=X*0; for r=1:numel(Jr); Y=Y+Jr{r}*X*Jr{r}.'; end
end
function prog = jump_diag(prog,D0,D,DT,Jr,sl,dwl,gm,th,deg,n)
    if strcmp(dwl,'range')
        prog=mat_ge(prog,-(Jstar_ms(D0,Jr)-D)-sl*eye(n),gm,th,deg,n);   % all theta
    else
        prog=sosmatrixineq(prog,-(Jstar_ms(D0,Jr)-DT)-sl*eye(n));       % theta=Tmin
    end
end
function prog = jump_spn(prog,P0,P,PT,Jr,sys,g,dwl,gm,th,deg,n)
    Cd=mvst_getf(sys,'Cd',mvst_getf(sys,'Cc',mvst_getf(sys,'C',zeros(0,n))));
    Ed=mvst_getf(sys,'Ed',zeros(n,0)); Fd=mvst_getf(sys,'Fd',zeros(size(Cd,1),size(Ed,2)));
    md=size(Ed,2);
    if strcmp(dwl,'range')
        prog=spn_ge(prog,-(Jstar_ms(P0,Jr)-P+Cd.'*Cd),gm,th,deg,n);
        if md>0; prog=spn_ge(prog,g*eye(md)-(Ed.'*P0*Ed+Fd.'*Fd),gm,th,deg,md); end
    else
        prog=spn_pt(prog,-(Jstar_ms(P0,Jr)-PT+Cd.'*Cd),th,n);
        if md>0; prog=spn_pt(prog,g*eye(md)-(Ed.'*P0*Ed+Fd.'*Fd),th,md); end
    end
end

% ---- SOS / cone-membership helpers ------------------------------------ %
function prog = mat_ge(prog,M,gmul,th,deg,sz)
% M(tau) matrix-SOS (PSD) on the interval {gmul>=0}:  M = S0 + gmul*S1,
% S0,S1 matrix-SOS.  gmul is degree 2, so S1 carries degree deg-2.
    [prog,S]=sospolymatrixvar(prog,monomials(th,0:max(0,deg-2)),[sz sz],'symmetric');
    prog=sosmatrixineq(prog,M-gmul*S); prog=sosmatrixineq(prog,S);
end
function prog = spn_ge(prog,M,gmul,th,deg,sz)
% M(tau) in SPN on the interval gmul>=0:  M = M1 + M2, M1 matrix-SOS,
% M2 symmetric with every entry a nonnegative polynomial there.
    [prog,M1]=sospolymatrixvar(prog,monomials(th,0:deg),[sz sz],'symmetric');
    [prog,M2]=sospolymatrixvar(prog,monomials(th,0:deg),[sz sz],'symmetric');
    prog=mat_ge(prog,M1,gmul,th,deg,sz);                     % M1 PSD on interval
    for i=1:sz; for j=i:sz
        [prog,s1]=sossosvar(prog,monomials(th,0:deg));
        prog=sosineq(prog, M2(i,j)-gmul*s1);                 % M2_ij >= 0 on interval
    end; end
    prog=soseq(prog, M-M1-M2);                               % M = M1+M2
end
function prog = entrywise_ge(prog,M,gmul,th,deg,sz)
% Every entry of the symmetric M(tau) is a nonnegative polynomial on the
% interval {gmul>=0}:  M_ij = s0 + gmul*s1, s0,s1 scalar SOS.
    for i=1:sz; for j=i:sz
        [prog,s1]=sossosvar(prog,monomials(th,0:deg));
        prog=sosineq(prog, M(i,j)-gmul*s1);
    end; end
end
function prog = spn_pt(prog,M,th,sz)
% M in SPN at a fixed timer value (constant / decision-affine matrix):
% M = M1 + M2, M1 PSD, M2 symmetric entrywise >= 0 (def:pos:cones).
    [prog,M1]=sospolymatrixvar(prog,monomials(th,0),[sz sz],'symmetric');
    [prog,M2]=sospolymatrixvar(prog,monomials(th,0),[sz sz],'symmetric');
    prog=sosmatrixineq(prog,M1);
    for i=1:sz; for j=i:sz; prog=sosineq(prog,M2(i,j)); end; end
    prog=soseq(prog,M-M1-M2);
end
function D = mkdiag(v,n)
% diagonal polynomial matrix with the polynomial vector v on the diagonal.
    D = v(1)*eeij(n,1);
    for i=2:n; D = D + v(i)*eeij(n,i); end
end
function E = eeij(n,i); E=zeros(n); E(i,i)=1; end
function c = tocell(X)
    if iscell(X); c=X; elseif isempty(X); c={}; else; p=size(X,3); c=cell(1,p);
        for k=1:p; c{k}=X(:,:,k); end; end
end
