function md = mvst_tothin(sys, varargin)
%MVST_TOTHIN  Convert a rich SOS-layer struct back to a thin mode struct.
%   MD = MVST_TOTHIN(RICH) maps a single-mode rich struct (SOS naming, drift
%   .A0, input .B0, noise .Ac/.Bc, jump/reset .J/.BJ/.Jf/.lamf/.BJf) back to
%   the thin OC/Riccati convention of MVST_MODE (drift .A, input .B, ...).
%   It is the inverse of MVST_TORICH on the fields both conventions share.
%   See lib/MODEL_CONVENTION.md.
%
%   MD = MVST_TOTHIN(RICH,'Q',Q,'R',R,'S',S,'W',W,'X0',X0,'ZT',ZT) supplies
%   the quadratic-cost data, which the rich struct DOES NOT CARRY; anything
%   omitted takes the MVST_MODE default (Q=0, R=I, S=0, W=0, X0=I).
%
%   A struct that is already thin (has .A, no .A0) is returned unchanged
%   (idempotent).  Fields that have no thin counterpart make the conversion
%   ill-posed and raise an error listing them:
%     * mode-indexed data (.A0 with M>1 modes, or a 4-D .Ac): a thin mode is
%       a SINGLE mode, so assemble the modes separately and wrap them with
%       MVST_SYSTEM instead;
%     * dissipativity input/output blocks (.Bw,.Dwc,.C,.Dw,.Bdf,.Bd,.CJ,.DJ)
%       and the generator/transition couplings (.Pig,.P): these describe the
%       supply-rate / coupling structure the thin mode cannot express.
%
%   See also MVST_TORICH, MVST_MODE, MVST_SYSTEM.

if ~isstruct(sys)
    error('mvst:tothin:input','SYS must be a struct.');
end

% Already thin: nothing to do (idempotent).
if isfield(sys,'A') && ~isfield(sys,'A0')
    md = sys; return;
end
if ~isfield(sys,'A0')
    error('mvst:tothin:input','SYS is not a rich struct (no .A0 field).');
end

% ---- collect unrepresentable fields ------------------------------------
bad = {};
if size(sys.A0,3) > 1;  bad{end+1} = 'A0 (mode-indexed, M>1)'; end
if isfield(sys,'Ac') && ndims(sys.Ac) > 3; bad{end+1} = 'Ac (mode-indexed, 4-D)'; end
noHome = {'Bw','Dwc','C','Dw','Bdf','Bd','CJ','DJ','Pig','P'};
for f = noHome
    if isfield(sys,f{1}) && ~isempty(sys.(f{1})); bad{end+1} = f{1}; end %#ok<AGROW>
end
if ~isempty(bad)
    error('mvst:tothin:unrepresentable', ...
        ['Cannot convert to a thin mode; these rich fields have no thin ', ...
         'counterpart:\n    %s\n(mode-indexed data: assemble each mode and ', ...
         'wrap with MVST_SYSTEM; dissipativity I/O: keep the rich struct ', ...
         'for the SOS layer.)'], strjoin(bad, ', '));
end

% ---- cost data (not carried by the rich struct) ------------------------
n = size(sys.A0,1);
p = inputParser;
p.addParameter('Q',  []); p.addParameter('R',  []); p.addParameter('S',  []);
p.addParameter('W',  []); p.addParameter('X0', []); p.addParameter('ZT', []);
p.parse(varargin{:});
o = p.Results;

A = sys.A0;
if isfield(sys,'B0') && ~isempty(sys.B0); B = sys.B0; else; B = zeros(n,0); end
m = size(B,2);
if isempty(o.Q); Q = zeros(n); else; Q = o.Q; end
if isempty(o.R); R = eye(m);   else; R = o.R; end

% ---- build the thin mode, forwarding the jump/reset extensions ----------
args = {};
if isfield(sys,'Ac')  && ~isempty(sys.Ac);   args = [args, {'Ac',  sys.Ac}];  end
if isfield(sys,'Bc')  && ~isempty(sys.Bc);   args = [args, {'Bc',  sys.Bc}];  end
if ~isempty(o.S);  args = [args, {'S',  o.S}];  end
if ~isempty(o.W);  args = [args, {'W',  o.W}];  end
if ~isempty(o.X0); args = [args, {'X0', o.X0}]; end
if isfield(sys,'J')   && ~isempty(sys.J);     args = [args, {'J',   sys.J}];   end
if isfield(sys,'BJ')  && ~isempty(sys.BJ);    args = [args, {'BJ',  sys.BJ}];  end
if isfield(sys,'Jf')  && ~isempty(sys.Jf);    args = [args, {'Jf',  sys.Jf}];  end
if isfield(sys,'lamf')&& ~isempty(sys.lamf);  args = [args, {'lamf',sys.lamf}];end
if isfield(sys,'BJf') && ~isempty(sys.BJf);   args = [args, {'BJf', sys.BJf}]; end

md = mvst_mode(A, B, Q, R, args{:});
if ~isempty(o.ZT); md.ZT = o.ZT; end
end
