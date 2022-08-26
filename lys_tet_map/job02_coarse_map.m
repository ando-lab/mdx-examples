%% Job 2: generate coarse maps

useParallel = true; % run in parallel (requires more RAM and parallel computing toolbox)

%%
% Filter, Integrate, and Correct

workingDir = fullfile('proc',{'2','3','4','5','6','7','8','9'});
%%
opts = struct(...
    'workingDirectory',workingDir,...
    'ndiv',[3,3,5],...  % this choice impacts speed to a large extent
    'window',2,...
    'maxCount',20,...
    'smax',Inf,...
    'binMode','coarse',...
    'minimumCounts',10,...
    'minimumPixels',10,...
    'parallel',useParallel,...
    'binExcludedVoxels',true...  % treats excluded voxels as Bragg peaks
    );

% assign run options
for j=1:length(opts)
    opts(j).run = {'filter','integrate','correct','export'};
end

[tf,EM] = proc.Batch.autorun(opts);

%% combine wedges

exportIn = fullfile({'2','3','4','5','6','7','8','9'},'export.mat');

[tf,EM] = proc.Batch.combine(...
    'workingDirectory','proc',...
    'exportIn',exportIn,...
    'combineBragg',true,...
    'mergeNeighbors',true);
