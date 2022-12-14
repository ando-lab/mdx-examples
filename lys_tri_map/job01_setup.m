%%

%% make directories

mkdir proc
mkdir export
mkdir images

%% Fetch diffraction images from SBGrid databank

!rsync -av rsync://data.sbgrid.org/10.15785/SBGRID/747/ ./images/

%% Check that transfer was successful

!cd images ; shasum -c files.sha

%% Set paths for importing geometry from XDS and image headers

% Create processing directories
mkdir proc

workingDir = fullfile('proc',...
    {'7_1','8_1','8_2','8_3','9_1','9_2',...
    '10_1','10_2','10_3','10_6','10_7'});

for j=1:numel(workingDir)
    mkdir(workingDir{j});
end

% Set relative paths for imports

xdsDir = fullfile('../../import',...
    {'xds_7_1','xds_8_1','xds_8_2', 'xds_8_3','xds_9_1','xds_9_2',...
    'xds_10_1','xds_10_2','xds_10_3','xds_10_6','xds_10_7'});

% Background frame information

backgroundTemplates = fullfile('images',...
        {'lys_nitr_7_bkg_1_????.cbf',...
        'lys_nitr_8_bkg_1_????.cbf',...
        'lys_nitr_8_bkg_1_????.cbf',...
        'lys_nitr_8_bkg_1_????.cbf',...
        'lys_nitr_9_bkg_2_????.cbf',...
        'lys_nitr_9_bkg_2_????.cbf',...
        'lys_nitr_10_bkg_1_????.cbf',...
        'lys_nitr_10_bkg_1_????.cbf',...
        'lys_nitr_10_bkg_1_????.cbf',...
        'lys_nitr_10_bkg_1_????.cbf',...
        'lys_nitr_10_bkg_1_????.cbf'});

backgroundFrameRanges = {...
        [1,50],...
        [171,220],...
        [126,175],...
        [1,50],...
        [1,50],...
        [46,95],...
        [251,300],...
        [296,345],...
        [1,360],...
        [71,120],...
        [116,165]};

%% import geometry

opts = struct(...
    'workingDirectory',workingDir,...
    'xdsDir',xdsDir,...
    'fileNameTemplate',backgroundTemplates,...
    'frameRange',backgroundFrameRanges);

% assign run options
for j=1:length(opts)
    opts(j).run = {'xds2geom','cbf2geom'};
end

[tf,EM] = proc.Batch.autorun(opts);

clear opts

%% define unit cell contents for scaling

%%
% *define unit cell contents*

[~,~,C] = io.mtz.read('import/XDS_ASCII_aimless.mtz');

SpaceGroup = symm.SpaceGroup(C.spaceGroupNumber);
asuPerUC = numel(SpaceGroup.generalPositions);
Crystal = geom.Crystal(C);
Vc = Crystal.UnitCell.vCell;
%%
lys_sequence = [...
    'KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINS',...
    'RWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDV',...
    'QAWIRGCRL'];
P = model.chem.Peptide(lys_sequence);
Lys = P.build();
%%
HOH = model.chem.Water(); % hard-coded (could use HOH.cif instead...)

data_NO3 = io.mmcif.searchFile('import/NO3.cif','NO3');
NO3 = io.mmcif.convert2Molecule(io.mmcif.parseRecord(data_NO3.record));

% *totals*

nHOH = 290*asuPerUC; % estimated
nNO3 = 6*asuPerUC;
nLys = 1*asuPerUC;

Molecules = [Lys,NO3,HOH];
occupancies = [nLys,nNO3,nHOH];

save('proc/unitCellInventory.mat',...
    'Molecules','occupancies','Crystal','Vc','asuPerUC');
