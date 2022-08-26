%% Job 2: generate coarse maps

useParallel = true; % run in parallel (requires more RAM and parallel computing toolbox)
ndiv = [5,5,11];

%%
% *generate fine grids for reintegration*


workingDir = fullfile('proc',{'2','3','4','5','6','7','8','9'});

opts = struct(...
    'workingDirectory',workingDir,...
    'matOut','grid.mat',...
    'smax',Inf,...
    'getSymmetryEquivalents',true,...
    'ndiv',ndiv,...
    'excludeBraggPosition',true);

for j=1:length(opts)
    [tf,EM] = proc.Batch.grid(opts(j));
end
%%
% *Re-integrate using fine grids*

opts = struct(...
    'workingDirectory',workingDir,...
    'geometryIn','geom.mat',...
    'bkgGeometryIn','geomBkg.mat',...
    'gridIn','grid.mat',...
    'scaleIn','scale.mat',...
    'matOut','reintegrate.mat',...
    'logOut','reintegrate.log',...
    'minimumCounts',0,...
    'minimumPixels',10,...
    'parallel',useParallel);

for j=1:length(opts)
    [tf,EM] = proc.Batch.reintegrate(opts(j));
end
%%
% *merge the result*

% get crystal info
load proc/unitCellInventory.mat Crystal

M = proc.script.MergeScaledDiffuse(...
    'Grid',grid.Sub3d('ndiv',ndiv),...
    'Crystal',Crystal,...
    'nMin',2,...
    'workingDirectory','proc');

files = fullfile({'2','3','4','5','6','7','8','9'},'reintegrate.mat');

fn = M.mapToColumns(files);

[hklTable,isincl] = M.mergeColumns(fn);

% split and merge for calculating cc1/2, cc*
rng(0,'twister'); % for reproducibility
[hklTable2] = M.mergeRandomHalfSets(fn,isincl);

M.clearTmp(fn);

% save table to mat file
hklMerge = hklTable;
save('proc/mergeFine.mat','hklMerge');
clear hklMerge hklTable isincl M

hklMerge = hklTable2;
save('proc/mergeFineSplit.mat','hklMerge');

clear hklMerge hklTable2
