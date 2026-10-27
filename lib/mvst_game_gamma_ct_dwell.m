function [gam, Kc, info] = mvst_game_gamma_ct_dwell(Ac, Bu, Bw, Q, Ru, Td, Chan, Dwc, Tlong, rtol, threl, ctol)
%MVST_GAME_GAMMA_CT_DWELL  Certified achievable H-infinity level, switched.
%   [GAM,KC] = MVST_GAME_GAMMA_CT_DWELL(AC,BU,BW,Q,RU,TD) returns the
%   smallest gamma for which (a) the per-mode game Riccati equations
%   admit central gains K_i at level gamma, and (b) the CLOSED-LOOP
%   switched system (drifts A_i + Bu K_i, performance output
%   z = [Q^{1/2} x; Ru^{1/2} u]) carries the MDT dissipativity
%   certificate of MVST_DISS_DWELL_CERT at the same level over every
%   schedule with dwell time at least TD.  GAM is thus a certified UPPER
%   BOUND on the achievable switched H-infinity level with mode-wise
%   central gains; the worst-mode critical level is a lower bound.  KC
%   returns the certified gains.
%
%   [...] = MVST_GAME_GAMMA_CT_DWELL(...,CHAN,DWC) carries the per-mode Brownian
%   diffusion channels CHAN (1 x M cell) and shared disturbance-multiplicative
%   channels DWC (1 x q cell) through the worst-mode lower bound and the per-mode
%   central-gain game Riccati (MVST_GAME_GAMMA_CT / MVST_GAME_DRE_CT), and folds
%   the closed-loop diffusion into the MDT certificate (MVST_DISS_DWELL_CERT).
%   With CHAN and DWC absent or empty the noise-free bound is returned unchanged.
%
%   [GAM,KC,INFO] = MVST_GAME_GAMMA_CT_DWELL(...,CHAN,DWC,TLONG,RTOL,THREL,CTOL)
%   switches to the LONG-HORIZON criterion of MVST_GAME_GAMMA_CT ( omit RTOL or pass it empty to keep the
%   former computation unchanged): the worst-mode lower bound is computed by
%   MVST_GAME_GAMMA_CT(...,TLONG,CHAN{i},DWC,RTOL,THREL), and at a trial level
%   the per-mode central gains are accepted only if each mode's game Riccati,
%   integrated backward over [0,TLONG] from zero, stays finite with
%   ||P_i|| < 1e5 AND keeps the uniform Isaacs margin
%   lambda_min(GAMMA^2 I - sum_l DWC{l}'P_i DWC{l}) >= THREL*GAMMA^2 on the
%   whole horizon.  CTOL is the relative bracket of the certificate bisection
%   (default 5e-3, the former value).  INFO returns the worst-mode lower bound
%   and the final certificate bracket.
%
%   See also MVST_GAME_GAMMA_CT, MVST_DISS_DWELL_CERT, MVST_L2GAIN_CT_DWELL.
if ~iscell(Ac)
    error('mvst:game_gamma_ct_dwell:modes', 'AC must be a 1 x M cell of drifts.');
end
if size(Q,1) ~= size(Q,2)
    error('mvst:game_gamma_ct_dwell:notSquare', 'Q must be square.');
end
M = numel(Ac);
if nargin < 7; Chan = {}; end
if nargin < 8; Dwc  = {}; end
if isempty(Chan); chi = repmat({{}}, 1, M); else; chi = Chan; end
if nargin < 9  || isempty(Tlong); Tlong = 30; end
if nargin < 10; rtol = []; end
if nargin < 11 || isempty(threl); threl = 1e-6; end
if nargin < 12 || isempty(ctol); ctol = 5e-3; end
longMode = ~isempty(rtol);
Dwcc = Dwc; if ~iscell(Dwcc); Dwcc = {Dwcc}; end
sqQ = sqrtm(Q); sqR = sqrtm(Ru);
lo = 0;
for ii = 1:M
    if longMode
        lo = max(lo, mvst_game_gamma_ct(Ac{ii}, Bu, Bw, Q, Ru, Tlong, chi{ii}, Dwc, rtol, threl));
    else
        lo = max(lo, mvst_game_gamma_ct(Ac{ii}, Bu, Bw, Q, Ru, 30, chi{ii}, Dwc));
    end
end
loBound = lo;
hi = max(2*lo, 1.0); Kc = [];
for kk = 1:40
    [okh, Kh] = cert(hi);
    if okh; Kc = Kh; break; end
    hi = hi*2;
    if kk == 40
        warning('mvst:game_gamma_ct_dwell:nohi','No certifiable gamma bracketed.');
        gam = hi; return;
    end
end
lo = max(lo, 1e-6);
while (hi - lo) > ctol*max(1,hi)
    mid = 0.5*(lo+hi);
    [okm, Km] = cert(mid);
    if okm; hi = mid; Kc = Km; else; lo = mid; end
end
gam = hi;
info = struct('lower',loBound,'lo',lo,'hi',hi,'ctol',ctol,'long',longMode);

    function [ok, K] = cert(g)
        K = cell(1,M); Acl = cell(1,M); Ccl = cell(1,M);
        for i = 1:M
            try
                [tt, Pt, Kut] = mvst_game_dre_ct(Ac{i}, Bu, Bw, Q, Ru, g, [0 Tlong], ...
                                                 [], [], chi{i}, Dwc);
            catch
                ok = false; return
            end
            if norm(Pt(:,:,1)) > 1e5; ok = false; return; end
            if longMode
                % boundedness over the whole horizon and uniform Isaacs margin
                if tt(1) > 1e-6 || ~all(isfinite(Pt(:))); ok = false; return; end
                for k = 1:size(Pt,3)
                    Th = g^2*eye(size(Bw,2));
                    for l = 1:numel(Dwcc); Th = Th - Dwcc{l}.'*Pt(:,:,k)*Dwcc{l}; end
                    if norm(Pt(:,:,k)) > 1e5 || min(eig((Th+Th.')/2)) < threl*g^2
                        ok = false; return
                    end
                end
            end
            K{i} = Kut(:,:,1);
            Acl{i} = Ac{i} + Bu*K{i};
            Ccl{i} = [sqQ; sqR*K{i}];
        end
        Dcl = zeros(size(Ccl{1},1), size(Bw,2));
        % The closed-loop state diffusion is unchanged by state feedback (no
        % input-multiplicative channel in the dwell game); pass CHI / DWC through.
        ok = mvst_diss_dwell_cert(Acl, Bw, Ccl, Dcl, Td, g, chi, Dwc);
    end
end
