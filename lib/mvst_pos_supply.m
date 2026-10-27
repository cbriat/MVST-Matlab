function res = mvst_pos_supply(sys, opts)
%MVST_POS_SUPPLY  General-linear-supply dissipativity by copositive linear
%   programming for LTI positive systems (Chapter "Linear Positive Systems"
%   of the monograph).  Tests dissipativity of the positive system with
%   respect to the general linear supply
%       sigma_f = q_c' z_c + r_c' w_c,   sigma_d = q_d' z_d + r_d' w_d,
%   and, when dissipative, returns the linear storage V = lambda' x.
%
%   RES = MVST_POS_SUPPLY(SYS, OPTS) implements the copositive program
%   eq:pos:genlinLP of Prop. prop:pos:genlinsupply, in the LTI constant-data
%   case with a constant storage weight lambda >= 0:
%     flow:  A' lambda - C_c' q_c            <= 0,
%            E_c' lambda - F_c' q_c - r_c    <= 0,
%     jump:  (J' - I) lambda - C_d' q_d      <= 0,     (when a jump channel is
%            E_d' lambda - F_d' q_d - r_d    <= 0.      given)
%   The weights q_c, r_c, q_d, r_d are of arbitrary sign.  Setting
%   q = -1, r = gamma*1 recovers the trace/L1 gain program; setting
%   q = -1, r = 1 recovers the mass-balance (passivity) supply of
%   Prop. prop:pos:disspassive.
%
%   SYS FIELDS (continuous channel; thin or rich naming accepted):
%     .A  or .A0    Metzler drift                (n x n)
%     .Ec or .E     disturbance input            (n x mc)
%     .Cc or .C     performance output           (pc x n)
%     .Fc or .F     feedthrough                  (pc x mc; default 0)
%   SYS FIELDS (jump channel; optional, all four enable the jump rows):
%     .J            reset map                     (n x n)
%     .Ed           jump disturbance input        (n x md)
%     .Cd           jump output                   (pd x n)
%     .Fd           jump feedthrough              (pd x md; default 0)
%
%   OPTS FIELDS:
%     .supply   struct with the supply weights (required):
%                 .qc (pc x 1), .rc (mc x 1), and, for the jump channel,
%                 .qd (pd x 1), .rd (md x 1).
%               Scalars are broadcast to the right length; a bare scalar in
%               .qc/.rc etc. means that value on every channel.
%     .margin   optional slack tightening the inequalities to <= -margin
%               (default 0, i.e. the book's non-strict dissipation form).
%     .lammin   lower bound on the storage weight (default 0).
%     plus any field consumed by MVST_POS_LP (.solver, .lpalgo, .lpopts).
%
%   RES FIELDS:
%     .feas     logical: the system is dissipative w.r.t. the supply.
%     .lambda   the storage weight (n x 1), so V = lambda' x; [] if infeas.
%     .lp       the LP data actually solved (LINPROG naming; see MVST_POS_LP).
%     .hasjump  logical: whether the jump rows were included.
%
%   Example (mass-balance passivity of a compartmental system):
%     sys = struct('A',[-0.9 0;0.5 -0.6],'Ec',[1;0],'Cc',[1 1]);
%     sup = struct('qc',-1,'rc',1);           % q = -1, r = 1
%     r = mvst_pos_supply(sys, struct('supply',sup));   % r.feas true if diss.
%
%   See also MVST_POS_STABILITY, MVST_POS_GAIN, MVST_POS_STABILIZE,
%   MVST_POS_LP.

    arguments
        sys  (1,1) struct
        opts (1,1) struct = struct()
    end

    sup = mvst_getf(opts, 'supply', []);
    if isempty(sup)
        error('mvst:pos:supply', ...
              'OPTS.supply must carry the weights .qc, .rc (and .qd, .rd).');
    end
    margin = mvst_getf(opts, 'margin', 0);
    lammin = mvst_getf(opts, 'lammin', 0);

    A  = mvst_getf(sys, 'A',  mvst_getf(sys, 'A0', []));
    Ec = mvst_getf(sys, 'Ec', mvst_getf(sys, 'E', []));
    Cc = mvst_getf(sys, 'Cc', mvst_getf(sys, 'C', []));
    if isempty(A) || isempty(Ec) || isempty(Cc)
        error('mvst:pos:supplydata', ...
              'SYS must carry the drift .A, disturbance .Ec and output .Cc.');
    end
    n  = size(A,1);
    mc = size(Ec,2);
    pc = size(Cc,1);
    Fc = mvst_getf(sys, 'Fc', mvst_getf(sys, 'F', zeros(pc, mc)));
    if isscalar(Fc) && (pc ~= 1 || mc ~= 1); Fc = Fc*ones(pc, mc); end

    qc = bcast(mvst_getf(sup, 'qc', []), pc, 'supply.qc');
    rc = bcast(mvst_getf(sup, 'rc', []), mc, 'supply.rc');

    % ---- flow rows ---------------------------------------------------- %
    Aub = [ A.'  ;      % A'lambda - Cc'qc <= 0
            Ec.' ];     % Ec'lambda - Fc'qc - rc <= 0
    bub = [ Cc.'*qc ;
            Fc.'*qc + rc ];

    % ---- jump rows (optional) ----------------------------------------- %
    J  = mvst_getf(sys, 'J',  []);
    Cd = mvst_getf(sys, 'Cd', []);
    Ed = mvst_getf(sys, 'Ed', []);
    hasjump = ~isempty(J) && ~isempty(Cd) && ~isempty(Ed);
    if hasjump
        md = size(Ed,2);
        pd = size(Cd,1);
        Fd = mvst_getf(sys, 'Fd', zeros(pd, md));
        if isscalar(Fd) && (pd ~= 1 || md ~= 1); Fd = Fd*ones(pd, md); end
        qd = bcast(mvst_getf(sup, 'qd', []), pd, 'supply.qd');
        rd = bcast(mvst_getf(sup, 'rd', []), md, 'supply.rd');
        Aub = [ Aub ;
                (J.' - eye(n)) ;   % (J'-I)lambda - Cd'qd <= 0
                Ed.' ];            % Ed'lambda - Fd'qd - rd <= 0
        bub = [ bub ;
                Cd.'*qd ;
                Fd.'*qd + rd ];
    end

    bub = bub - margin;            % optional strictening

    lp = struct('nvar', n, 'Aub', Aub, 'bub', bub, 'lb', lammin*ones(n,1));
    o  = mvst_pos_lp(lp, opts);

    res = struct('feas', o.feas, 'lambda', [], 'lp', o.lp, 'hasjump', hasjump);
    if o.feas
        res.lambda = o.x;
    end
end

% --------------------------------------------------------------------- %
function v = bcast(v, k, name)
    if isempty(v)
        error('mvst:pos:supply', 'OPTS.%s is required.', name);
    end
    if isscalar(v); v = v*ones(k,1); return; end
    v = v(:);
    if numel(v) ~= k
        error('mvst:pos:supply', ...
              'OPTS.%s must be scalar or length %d (got %d).', name, k, numel(v));
    end
end
