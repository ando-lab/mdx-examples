%%

%% make directories

mkdir proc
mkdir export
mkdir images

%% Fetch diffraction images from SBGrid databank

!rsync -av rsync://data.sbgrid.org/10.15785/SBGRID/958/ ./images/

%% Check that transfer was successful

!cd images ; shasum -c files.sha

%% Set paths for importing geometry from XDS and image headers

% Create processing directories
mkdir proc

workingDir = fullfile('proc',...
    {'2_38','2_39','2_40','2_41','5_51','5_52',...
    '5_53','5_54'});

for j=1:numel(workingDir)
    mkdir(workingDir{j});
end

% Set relative paths for imports

xdsDir = fullfile('../../import',...
    {'xds_2_38','xds_2_39','xds_2_40','xds_2_41','xds_5_51','xds_5_52',...
    'xds_5_53','xds_5_54'});

% Background frame information

backgroundTemplates = fullfile('images',...
        {'lys_rt_2_42_????.cbf',...
        'lys_rt_2_42_????.cbf',...
        'lys_rt_2_42_????.cbf',...
        'lys_rt_2_42_????.cbf',...
        'lys_rt_5_55_????.cbf',...
        'lys_rt_5_55_????.cbf',...
        'lys_rt_5_55_????.cbf',...
        'lys_rt_5_55_????.cbf'});

backgroundFrameRanges = {...
        [1,50],...
        [91,140],...
        [181,230],...
        [271,320],...
        [1,50],...
        [91,140],...
        [181,230],...
        [271,320]};

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

data_CL = io.mmcif.searchFile('import/Cl.cif','CL');
Cl = io.mmcif.convert2Molecule(io.mmcif.parseRecord(data_CL.record));
% *totals*

nHOH = 472.5*asuPerUC; % estimated from MD simulation
nCl = 8*asuPerUC;
nLys = 1*asuPerUC;

Molecules = [Lys,Cl,HOH];
occupancies = [nLys,nCl,nHOH];

save('proc/unitCellInventory.mat',...
    'Molecules','occupancies','Crystal','Vc','asuPerUC');
