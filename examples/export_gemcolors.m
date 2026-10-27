%% EXPORT_GEMCOLORS  Write the gem palette to a LaTeX \definecolor file.
%  Produces figures/data/gemcolors.tex so the TikZ diagrams use exactly the
%  same colors as the generated Matlab figures.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
M = mvst_gemcolors().order;
d = fullfile(mvst_figdir(),'data');
if ~exist(d,'dir'); mkdir(d); end
fid = fopen(fullfile(d,'gemcolors.tex'),'w');
names = {'gemblue','gemred','gemyellow','gempurple','gemgreen','gemcyan','gemmaroon'};
for i = 1:min(numel(names), size(M,1))
    fprintf(fid, '\\definecolor{%s}{rgb}{%.4f,%.4f,%.4f}\n', names{i}, M(i,1), M(i,2), M(i,3));
end
fclose(fid);
fprintf('wrote %s\n', fullfile(d,'gemcolors.tex'));
