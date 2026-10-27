function [nu, P] = mvst_isp_margin_ct(A, Bw, C, Dw, Cc, Dwc, opts)
%MVST_ISP_MARGIN_CT  Input-strict-passivity index of a (stochastic) LTI flow.
%   NU = MVST_ISP_MARGIN_CT(A,BW,C,DW,CC,DWC) returns the largest index nu of
%   input strict passivity, in the convention of the Dissipativity chapter
%   (Section sec:diss:passivity: supply He(w'z) - nu |w|^2, i.e. the passivity
%   supply eq:diss:M:pass with M22 <- He(Dw) - nu I), for the Ito system
%       dx = (A x + Bw w)dt + sum_l (Cc{l} x + Dwc{l} w) dW_l,   z = C x + Dw w.
%   CC (state-multiplicative) and DWC (disturbance-multiplicative) are 1 x q
%   cells (either may be empty; the shorter is padded with zeros).
%
%   NU is the optimal value of the SDP, in the constant storage matrix P
%   (P = -S, S >= 0 the storage) and nu,
%       maximize nu  subject to  -P >= 0  and
%       [ A'P + P A + sum Cc'P Cc ,  P Bw + sum Cc'P Dwc + C'           ]
%       [        (.)'             ,  He(Dw) - nu I + sum Dwc'P Dwc      ] >= 0,
%   the constant-storage flow condition C_M(P) >= 0 of the chapter for the
%   ISP supply.  It is solved with SOSTOOLS/SeDuMi (constant certificate, so
%   the SOS program is the LMI itself).  Without noise, NU equals the
%   frequency-domain value min_w lambda_min(G(jw)+G(jw)'), the same convention
%   (G+G' >= nu I) as MVST_PASSIVITY_CT and MVST_PASSIVITY since 2026-10-27.
%
%   OPTS.solver selects the SDP backend (default 'sedumi').  P is returned.
%
%   See also MVST_PASSIVITY_CT, MVST_DISS_SOS.
if nargin < 5 || isempty(Cc);  Cc  = {}; end
if nargin < 6 || isempty(Dwc); Dwc = {}; end
if nargin < 7; opts = struct(); end
if ~iscell(Cc);  Cc  = {Cc};  end
if ~iscell(Dwc); Dwc = {Dwc}; end
n = size(A,1); m = size(Bw,2);
assert(size(C,1) == m, 'mvst:isp_margin_ct:square', ...
    'Passivity needs a square port (size(C,1) == size(Bw,2)).');
q = max(numel(Cc), numel(Dwc));
Cc(end+1:q)  = {zeros(n)};
Dwc(end+1:q) = {zeros(n,m)};
oo = mvst_opts(opts, true, 'sos'); solver = oo.solver;
Rcap = 1e4;                                   % storage bound 0 <= -P <= Rcap I

pvar th; nuv = dpvar('nuv');
prog = sosprogram(th); prog = sosdecvar(prog, nuv);
[prog, Pv] = sospolymatrixvar(prog, monomials(th,0), [n n], 'symmetric');
G11 = A.'*Pv + Pv*A; G12 = Pv*Bw + C.'; G22 = (Dw + Dw.') - nuv*eye(m);
for l = 1:q
    G11 = G11 + Cc{l}.'*Pv*Cc{l};
    G12 = G12 + Cc{l}.'*Pv*Dwc{l};
    G22 = G22 + Dwc{l}.'*Pv*Dwc{l};
end
CM = [G11, G12; G12.', G22];
prog = sosmatrixineq(prog, CM);
prog = sosmatrixineq(prog, -Pv);
prog = sosmatrixineq(prog, Rcap*eye(n) + Pv);
prog = sossetobj(prog, -nuv);
o = mvst_sosopts(oo);
prog = mvst_sossolve(prog, o);
nu = double(sosgetsol(prog, nuv));
P  = double(sosgetsol(prog, Pv));
end
