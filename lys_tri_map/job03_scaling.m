%% Job 3: refine scaling model to coarse maps

% add mdx-lib to path
addpath('~/Documents/GitHub/ando-lab/mdx-lib')

%%

[tf,EM] = proc.Batch.scale(...
    'workingDirectory','proc',...
    'scaleIn','combine.mat',...
    'csMult',10,'cizMult',1); % relaxing regularization of c improved cc 1/2

% *merge Diffuse intensity*
[tf,EM] = proc.Batch.merge(...
    'workingDirectory','proc');

% *merge Bragg intensity*
[tf,EM] = proc.Batch.merge(...
    'workingDirectory','proc',...
    'matOut','mergeBragg.mat',...
    'logOut','mergeBragg.log',...
    'mergeBragg',true,...
    'sigmaCutoff',2); % be harsh with outliers!

%%
% *replace Bragg intensities by those integrated with xds / aimless*

mergeBragg = load('proc/mergeBragg.mat','hklMerge');

S2R = proc.script.ScaleToReference(...
    'mtzIn','import/XDS_ASCII_aimless.mtz',...
    'refcols',{'hasu'  'kasu'  'lasu'  'I'  'sigma'},...
    'outcols',{'hasu'  'kasu'  'lasu'  'I'  'sigma'},... % use same column names in output table
    'nIter',1E4); % maximum number of iterations for outlier rejection.

% load aimless and scale 2 bragg
[hklMerge,hklRef,scaleFactor] = S2R.run(mergeBragg.hklMerge);

save('proc/mergeBraggAimless.mat','hklMerge','hklRef','scaleFactor');

clear hklMerge hklRef scaleFactor S2R mergeBragg

%%
% *absolute intensity scaling by modofied Krogh-Moe method*

[tf,EM] = proc.Batch.rescale(...
    'workingDirectory','proc/',...
    'unitCellInventoryIn','unitCellInventory.mat',...
    'mergeIn','merge.mat',... 
    'mergeBraggIn','mergeBraggAimless.mat',... 
    'mergeOut','mergeAbsolute.mat',... % bragg and diffuse together
    'scaleIn','scale.mat',...
    'scaleOut','scaleAbsolute.mat',...
    'matOut','rescale.mat',... % stats about scaling are saved here
    'logOut','rescale.log',...
    'smax',1,...
    'npts',501,...
    'scutoff',.7); % maximum 1/d to use in scaling

%%
% *copy scaling models to wedge directories*

% get list of wedge scatterbrain directories from combine.mat
load('proc/combine.mat','options');
exportDir = cellfun(@fileparts,options.exportIn,'Uni',0);

% get scaling model
scale = load('proc/scaleAbsolute.mat');

% export scaling models
for j=1:length(exportDir)
    ScalingModel = scale.ScalingModel(j);
    save(fullfile('proc/',exportDir{j},'scale.mat'),'ScalingModel');
end