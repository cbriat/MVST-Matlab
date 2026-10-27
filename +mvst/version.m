function v = version()
%MVST.VERSION  Version string of the Matrix-Valued Systems Toolbox (MVST).
%   V = MVST.VERSION returns the toolbox version as a character vector,
%   for example '1.0.0-rc1' (semantic versioning; see CHANGELOG.md).
%
%   Example:
%       fprintf('MVST %s\n', mvst.version);
%
%   See also MVST.SOLVERCHECK, MVST.SYSTEM.
    v = '1.0.0-rc1';
end
