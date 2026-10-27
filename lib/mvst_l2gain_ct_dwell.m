function gam = mvst_l2gain_ct_dwell(Ac, Bw, Cc, Dw, Td, Chan, Dwc)
%MVST_L2GAIN_CT_DWELL  Dwell-time-certified L2-gain bound, switched system.
%   GAM = MVST_L2GAIN_CT_DWELL(AC,BW,CC,DW,TD) returns the smallest gamma
%   for which the MDT dissipativity certificate of MVST_DISS_DWELL_CERT
%   could be constructed: a certified UPPER BOUND on the L2-gain of the
%   switched system over every schedule with dwell time at least TD.  The
%   worst-mode gain is a lower bound; the true switched gain lies between
%   the two.
%
%   GAM = MVST_L2GAIN_CT_DWELL(...,CHAN,DWC) carries the per-mode Brownian
%   diffusion channels CHAN (1 x M cell) and shared disturbance-multiplicative
%   channels DWC (1 x q cell) through both the worst-mode lower bound
%   (MVST_L2GAIN_CT) and the MDT certificate (MVST_DISS_DWELL_CERT).  With
%   CHAN and DWC absent or empty the noise-free bound is returned unchanged.
%
%   See also MVST_DISS_DWELL_CERT, MVST_L2GAIN_CT, MVST_GAME_GAMMA_CT_DWELL.
if ~iscell(Ac) || ~iscell(Cc)
    error('mvst:l2gain_ct_dwell:modes', 'AC and CC must be 1 x M cells.');
end
if ~(isscalar(Td) && isreal(Td) && Td > 0)
    error('mvst:l2gain_ct_dwell:dwell', 'TD (dwell-time) must be a positive scalar.');
end
if nargin < 6; Chan = {}; end
if nargin < 7; Dwc  = {}; end
M = numel(Ac);
if isempty(Chan); chi = repmat({{}}, 1, M); else; chi = Chan; end
lo = 0;
for i = 1:M
    lo = max(lo, mvst_l2gain_ct(Ac{i}, Bw, Cc{i}, Dw, chi{i}, Dwc));   % worst mode
end
hi = max(2*lo, 1.0);
for kk = 1:40
    if mvst_diss_dwell_cert(Ac, Bw, Cc, Dw, Td, hi, Chan, Dwc); break; end
    hi = hi*2;
    if kk == 40
        warning('mvst:l2gain_ct_dwell:nohi','No certifiable gamma bracketed.');
        gam = hi; return;
    end
end
lo = max(lo, 1e-6);
while (hi - lo) > 5e-3*max(1,hi)
    mid = 0.5*(lo+hi);
    if mvst_diss_dwell_cert(Ac, Bw, Cc, Dw, Td, mid, Chan, Dwc)
        hi = mid;
    else
        lo = mid;
    end
end
gam = hi;
end
