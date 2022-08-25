%%

% add mdx-lib to path
addpath('~/Documents/GitHub/ando-lab/mdx-lib')

%% make directories

mkdir proc
mkdir export
mkdir images

%% Fetch diffraction images from SBGrid databank

!rsync -av rsync://data.sbgrid.org/10.15785/SBGRID/957/ ./images/

%% Check that transfer was successful

!cd images ; shasum -c files.sha

%% Set paths for importing geometry from XDS and image headers

% Create processing directories
mkdir proc

workingDir = fullfile('proc',...
    {'2','3','4','5','6','7',...
    '8','9'});

for j=1:numel(workingDir)
    mkdir(workingDir{j});
end

% Set relative paths for imports

xdsDir = fullfile('../../import',...
    {'xds_2','xds_3','xds_4', 'xds_5','xds_6','xds_7',...
    'xds_8','xds_9'});

% Background frame information

backgroundTemplates = fullfile('images','lys_1_bkg_1_????.cbf');
    
backgroundFrameRanges = [1,360];

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

lys_sequence = [...
    'KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINS',...
    'RWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDV',...
    'QAWIRGCRL'];
P = model.chem.Peptide(lys_sequence);
Lys = P.build();

HOH = model.chem.Water(); % hard-coded (could use HOH.cif instead...)

data_Cl = io.mmcif.searchFile('import/Cl.cif','CL');
Cl = io.mmcif.convert2Molecule(io.mmcif.parseRecord(data_Cl.record));

% *totals*

nHOH = 413*asuPerUC; % estimated from DAC simulation
nCl = 7*asuPerUC;
nLys = 1*asuPerUC;

Molecules = [Lys,Cl,HOH];
occupancies = [nLys,nCl,nHOH];

save('proc/unitCellInventory.mat',...
    'Molecules','occupancies','Crystal','Vc','asuPerUC');
