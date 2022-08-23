function Compute_Parameters_From_Lines

% List of the subjects
idFile     = '/data/hagmann_group/yaleman/OASIS_Data/allIds1.txt';

% Reding the IDs
Ids = textread(idFile,'%s');
Nsubj = length(Ids);

% addpath('/home/yalemang/matlab/scripts/ClusterScripts');
% addpath(genpath('/home/yalemang/matlab/Tools'));
% addpath('/home/yalemang/matlab/scripts/Surfaces_myTools');
% addpath('/home/yalemang/matlab/scripts/Miscellaneous_myTools');
%
% addpath('/home/yalemang/matlab/Tools/gptoolbox/mex');

fsdir = '/data/hagmann_group/yaleman/OASIS_Data/derivatives/freesurfer';
% function ApplyANTs_to_trkfiles

cont = 0;
for i = 1:Nsubj
    subjId = Ids{i};
    disp(['Processing subject: ' num2str(i)  ' of ' num2str(Nsubj)]);
    surfDir     = [fsdir filesep subjId filesep 'surf'];
    labelDir     = [fsdir filesep subjId filesep 'label'];
    
    lhpial  = [surfDir filesep 'lh.pial'];
    rhpial  = [surfDir filesep 'rh.pial'];
    slines  = [surfDir filesep 'sulcABLE.lines.morph.mat'];
    lhannot = [labelDir filesep 'lh.aparc.annot'];
    lhlobannot = [labelDir filesep 'lh.lobes+cing.annot'];
    rhannot = [labelDir filesep 'rh.aparc.annot'];
    rhlobannot = [labelDir filesep 'lh.lobes+cing.annot'];
    lhdestrieux = [labelDir filesep 'lh.aparc.a2009s.annot'];
    rhdestrieux = [labelDir filesep 'rh.aparc.a2009s.annot'];
    
    
    slength = [fsdir filesep subjId filesep 'stats' filesep subjId '_sulcABLE.lines.length.stats'];
    delete(slength);
    %% ========================= Left Hemisphere =========================== %%
    if exist(lhpial,'file')&exist(rhpial,'file')&exist(lhannot,'file')&exist(rhannot,'file')&exist(slines,'file')%&~exist(slength,'file')
        cont = cont +1;
        %% Left Hemisphere
        Pial = Read_Surface(lhpial);
        [~, boolparc, ctab] = Aparc2Lobes(lhannot, 1, lhlobannot);
        load(slines)
        sLinesIndex = morph.left.sulci;
        gLinesIndex = morph.left.gyri;
        
        [lhlabels,lhstructNumber,lhstructName] = Label_Sulci_Using_Aparc2009s(Pial,sLinesIndex,lhdestrieux);
        
        temp = boolparc(sLinesIndex);
        [a,b] = ismember(temp, ctab.table(:,5));
        
        Xp = [Pial.SurfData.vertices(sLinesIndex(:,1),1)'; Pial.SurfData.vertices(sLinesIndex(:,2),1)'];
        Yp = [Pial.SurfData.vertices(sLinesIndex(:,1),2)'; Pial.SurfData.vertices(sLinesIndex(:,2),2)'];
        Zp = [Pial.SurfData.vertices(sLinesIndex(:,1),3)'; Pial.SurfData.vertices(sLinesIndex(:,2),3)'];
        A = sqrt(sum(diff(Xp,1).^2 + diff(Yp,1).^2 + diff(Zp,1).^2,1));
        ind = find(b(:,1)~=1);
        lhVect = accumarray(b(ind,1),A(ind));
        lhVect = lhVect(2:end);
        
        %% Right Hemisphere
        Pial = Read_Surface(rhpial);
        [~, boolparc, ctab] = Aparc2Lobes(rhannot, 1, rhlobannot);
        sLinesIndex = morph.right.sulci;
        gLinesIndex = morph.right.gyri;
        
        temp = boolparc(sLinesIndex);
        [a,b] = ismember(temp, ctab.table(:,5));
        
        Xp = [Pial.SurfData.vertices(sLinesIndex(:,1),1)'; Pial.SurfData.vertices(sLinesIndex(:,2),1)'];
        Yp = [Pial.SurfData.vertices(sLinesIndex(:,1),2)'; Pial.SurfData.vertices(sLinesIndex(:,2),2)'];
        Zp = [Pial.SurfData.vertices(sLinesIndex(:,1),3)'; Pial.SurfData.vertices(sLinesIndex(:,2),3)'];
        A = sqrt(sum(diff(Xp,1).^2 + diff(Yp,1).^2 + diff(Zp,1).^2,1));
        ind = find(b(:,1)~=1);
        rhVect = accumarray(b(ind,1),A(ind));
        rhVect = rhVect(2:end);
        
        tempVar = [lhVect rhVect]';
        
        lengthVect = [sum(lhVect) + sum(rhVect);sum(lhVect);sum(rhVect);tempVar(:)];
        lhNames = strcat('lh-', ctab.struct_names(2:end));
        rhNames = strcat('rh-', ctab.struct_names(2:end));
        allNames = {'whole-brain';'lh-hemisphere';'rh-hemisphere'};
        tempNames = [lhNames rhNames]';
        allNames = [allNames;tempNames(:)];
        allNames = strcat( allNames,'-sulclength');
        resVar = [{'subjID', subjId };allNames num2cell(lengthVect)]';
        cad2textfile(cell2cad(resVar),slength);
        if cont==1
            allResults = resVar;
        else
            allResults = [allResults;resVar(2,:)];
        end
    end
end
cad2textfile(cell2cad(allResults),[fsdir filesep 'allsubj-sulcallength.csv']);

return;

function varargout = Aparc2Lobes(varargin);
%
% Syntax :
%     OutAnnotFile = Aparc2Lobes(AnnotFile, boolcingulate, OutAnnotFile);
%
% This function does an approximate mapping of individual 'Desikan-Killiany' ROIs
% (found in ?h.aparc.annot) to the lobes
%
% Input Parameters:
%        AnnotFile              : Annotation File
%        boolcingulate          : Boolean variable to include or not Cingulate Lobe
%                                 (0 = do not include, 1 = include)
%        OutAnnotFile           : Saving the lobar parcellation in an
%                                 annotation file
%
% Output Parameters:
%        OutAnnotFile           : Saving the lobar parcellation in an
%                                 annotation file. If nargin <3
%                                 OutAnnotFile is a vector file containing
%                                 lobar parcellation.
%
% See also: save_annotfiles
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% November 13th 2014
% Version $1.0

%% ============================= Checking Inputs ======================= %%
AnnotFile = varargin{1};
if ~exist(AnnotFile,'file')
    errordlg('The Annotation File Does not exist');
    return;
else
    try
        % Reading Parcellation
        [txt,colortable] = read_cfiles(AnnotFile);
        if ~isfield(colortable,'table')
            errordlg('The File is not an annotation file');
            return;
        end
    catch
        errordlg('Error reading Annotation File');
        return;
    end
end

if nargin == 1
    boolcingulate = 0;
end
if nargin == 2
    boolcingulate = logical(varargin{2});
end
if nargin == 3
    OutAnnotFile = varargin{3};
    boolcingulate = logical(varargin{2});
end
if nargin > 3
    errordlg('To Many Input Parameters');
    return;
end
% if nargout > 2
%     errordlg('To Many Output Parameters');
%     return;
% end
%% ========================= End of Checking Inputs ==================== %%

%% ============================= Defining Lobes ======================== %%
if boolcingulate
    % Frontal Lobe
    Fron = strvcat('caudalmiddlefrontal',...
        'lateralorbitofrontal',...
        'medialorbitofrontal',...
        'paracentral',...
        'parsopercularis',...
        'parsorbitalis',...
        'parstriangularis',...
        'precentral',...
        'rostralmiddlefrontal',...
        'superiorfrontal',...
        'frontalpole');
    
    % Parietal Lobe
    Par = strvcat('inferiorparietal',...
        'isthmuscingulate',...
        'postcentral',...
        'precuneus',...
        'superiorparietal',...
        'supramarginal');
    
    % Cingulate Lobe
    Cing = strvcat('rostralanteriorcingulate',...
        'caudalanteriorcingulate',...
        'posteriorcingulate',...
        'isthmuscingulate');
    
else
    % Frontal Lobe
    Fron = strvcat('caudalanteriorcingulate',...
        'caudalmiddlefrontal',...
        'lateralorbitofrontal',...
        'medialorbitofrontal',...
        'paracentral',...
        'parsopercularis',...
        'parsorbitalis',...
        'parstriangularis',...
        'precentral',...
        'rostralanteriorcingulate',...
        'rostralmiddlefrontal',...
        'superiorfrontal',...
        'frontalpole');
    
    % Parietal Lobe
    Par = strvcat('inferiorparietal',...
        'isthmuscingulate',...
        'postcentral',...
        'posteriorcingulate',...
        'precuneus',...
        'superiorparietal',...
        'supramarginal');
end

% Temporal Lobe
Temp = strvcat('bankssts',...
    'entorhinal',...
    'fusiform',...
    'inferiortemporal',...
    'middletemporal',...
    'parahippocampal',...
    'superiortemporal',...
    'temporalpole',...
    'transversetemporal');

% Occipital Lobe
Occ = strvcat('cuneus',...
    'lateraloccipital',...
    'lingual',...
    'pericalcarine');

% Insula Lobe
Ins = 'insula';
%% =======================End of Defining Lobes ======================== %%


%% ================= Removing white spaces from Structure Names ======== %%
totnames1 = char(colortable.struct_names);
totnames = '';
for i = 1:size(totnames1,1)
    totnames = strvcat(totnames,deblank(totnames1(i,:)));
end
ind = find(sum(isspace(totnames))==size(colortable.struct_names,1));
totnames(:,ind) = [];
%% ========== End of Removing white spaces from Structure Names ======== %%

%% ================= Removing white spaces from Frontal Lobe =========== %%
ind = find(sum(isspace(Fron))==size(Fron,1));
Fron(:,ind) = [];
ind = ismember(totnames(:,1:size(Fron,2)),Fron,'rows');
FronIds = colortable.table(ind,5);
indfron = ismember(txt,FronIds);
%% ========== End of Removing white spaces from Frontal Lobe =========== %%

%% ================= Removing white spaces from Parietal Lobe ========== %%
ind = find(sum(isspace(Par))==size(Par,1));
Par(:,ind) = [];
ind = ismember(totnames(:,1:size(Par,2)),Par,'rows');
ParIds = colortable.table(ind,5);
indpar = ismember(txt,ParIds);
%% ========== End of Removing white spaces from Parietal Lobe ========== %%

%% ================= Removing white spaces from Temporal Lobe ========== %%
ind = find(sum(isspace(Temp))==size(Temp,1));
Temp(:,ind) = [];
ind = ismember(totnames(:,1:size(Temp,2)),Temp,'rows');
TempIds = colortable.table(ind,5);
indtemp = ismember(txt,TempIds);
%% ========== End of Removing white spaces from Temporal Lobe ========== %%

%% ================= Removing white spaces from Occipital Lobe ========= %%
ind = find(sum(isspace(Occ))==size(Occ,1));
Occ(:,ind) = [];
ind = ismember(totnames(:,1:size(Occ,2)),Occ,'rows');
OccIds = colortable.table(ind,5);
indocc = ismember(txt,OccIds);
%% ========== End of Removing white spaces from Occipital Lobe ========= %%

if boolcingulate
    %% ================= Removing white spaces from Occipital Lobe ========= %%
    ind = find(sum(isspace(Cing))==size(Cing,1));
    Cing(:,ind) = [];
    ind = ismember(totnames(:,1:size(Cing,2)),Cing,'rows');
    CingIds = colortable.table(ind,5);
    indcing = ismember(txt,CingIds);
    %% ========== End of Removing white spaces from Occipital Lobe ========= %%
end

%% ============= Creating New Label Id vector and ColorTable =========== %%
lobparc = txt*0;

lobparc(indfron) = 1;  % Frontal Lobe
stnames{2,1} = 'frontallobe ';
lobparc(indpar) = 2;   % Parietal Lobe
stnames{3,1} = 'parietallobe ';
lobparc(indtemp) = 3;  % Temporal Lobe
stnames{4,1} = 'temporallobe ';
lobparc(indocc) = 4;   % Occipital Lobe
if boolcingulate
    lobparc(indcing) = 5;   % Cingulate Lobe
    stnames{5,1} = 'linguallobe ';
end

if find(unique(lobparc) == 0)
    stnames{1,1} = 'unknown';
    stnames{2,1} = 'frontallobe';
    lobparc(indpar) = 2;   % Parietal Lobe
    stnames{3,1} = 'parietallobe';
    lobparc(indtemp) = 3;  % Temporal Lobe
    stnames{4,1} = 'temporallobe';
    lobparc(indocc) = 4;   % Occipital Lobe
    stnames{5,1} = 'occipitallobe';
    if boolcingulate
        lobparc(indcing) = 5;   % Cingulate Lobe
        stnames{6,1} = 'linguallobe';
    end
else
    stnames{1,1} = 'frontallobe';
    lobparc(indpar) = 2;   % Parietal Lobe
    stnames{2,1} = 'parietallobe';
    lobparc(indtemp) = 3;  % Temporal Lobe
    stnames{3,1} = 'temporallobe';
    lobparc(indocc) = 4;   % Occipital Lobe
    stnames{4,1} = 'occipitallobe';
    if boolcingulate
        lobparc(indcing) = 5;   % Cingulate Lobe
        stnames{5,1} = 'linguallobe';
    end
end

sts = unique(lobparc);
if find(sts ==0)
    colors = [255 255 255; 12 73 153; 7 96 32; 210 113 26; 126 18 33; 187 205 31];
else
    colors = [12 73 153; 7 96 32; 210 113 26; 126 18 33; 187 205 31];
end
ctab = [colors colors(:,1)*0 colors(:,1)+colors(:,2)*2^8+colors(:,3)*2^16];
Nstruct = length(sts);
for i = 1:Nstruct
    ind = find(lobparc==sts(i));
    if sts(i) == 0
        ctab(i,:)  = [255 255 255 0 16777215];
        lobparc(ind) = 16777215;
    else
        lobparc(ind) = ones(size(ind,1),1)*ctab(i,5);
    end
end
colortable.numEntries = length(sts);
colortable.orig_tab = 'Custom Colortable';
colortable.struct_names = stnames;
colortable.table = ctab;

%% ====================== End of Defining Lobes ======================== %%

if nargin <= 2;
    [~, lobparc] = ismember(lobparc,colortable.table(:,5));
    varargout{1} = lobparc;
    varargout{2} = colortable;
    
else
    OutAnnotFile = save_annotfiles(lobparc,OutAnnotFile,colortable);
    varargout{1} = OutAnnotFile;
    varargout{2} = lobparc;
    varargout{3} = colortable;
end
return;
function OutAnnot = save_annotfiles(labelsd,OutAnnot,colortable);
%
% Syntax :
% OutAnnot = save_annotfiles(labelsd,OutAnnot,colortable);
%
% This script creates an annotation file for cortical parcelation.
%
%
% Input Parameters:
%   labelsd            : Structures Labels
%   OutAnnot           : New Annotation filename
%   colortable         : Freesurfer colortable in Matlab structure
%
% Output Parameters:
%   OutAnnot             : New Annotation filename
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% March 22th 2012
% Version $1.0

%=====================Checking Input Parameters===========================%
sts = unique(labelsd);
Nstruct = length(sts);
t = 1;
if nargin<3
    [colortable,labels] = Create_FS_Colortable(labelsd);
    ctab = colortable.table;
else
    labels = labelsd;
    ctab = colortable.table;
end
%=========================================================================%
%=====================     Main Program    ===============================%
Np = length(labels);
fid = fopen(OutAnnot,'wb','b');
fwrite(fid, Np, 'int');
Avalues(2:2:Np*2) = labels;
Avalues(1:2:Np*2) = [0:Np-1]';
fwrite(fid, Avalues, 'int');
fwrite(fid, 1, 'int');
fwrite(fid, -2, 'int');
fwrite(fid,length(unique(labels)),'int');
fwrite(fid, length(OutAnnot), 'int');
fwrite(fid, OutAnnot, '*char');
fwrite(fid,length(unique(labels)),'int');
for i = 1:Nstruct
    fwrite(fid, i-1, 'int');
    name = deblank(colortable.struct_names{i});
    len = length(name);
    fwrite(fid, len+1, 'int');
    fwrite(fid, [name ' '], 'char');
    fwrite(fid, ctab(i,1), 'int');
    fwrite(fid, ctab(i,2), 'int');
    fwrite(fid, ctab(i,3), 'int');
    fwrite(fid, ctab(i,4), 'int');
end
fclose(fid);
%=========================================================================%
return;

function varargout = cad2textfile(varargin);
%
% Syntax :
%   textFile = cad2textfile(cad2print, textFile);
%
% This script writes lines included from the variable cad2print into
% a text file.
%
% Input Parameters:
%         cad2print        : Lines to be printed.
%         textFile         : Output text file.
%
% Output Parameters:
%         functNames       :  Functions filenames
%
%
% Related references:
%
%
% See also: 
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% September 14th 2014
% Version $1.0

%% ============================= Checking Inputs ======================= %%


%% ========================== Checking Inputs ========================== %%
if nargin <2
    error('Two inputs are mandatory');
    return;
end
cad2print = varargin{1};
textFile = varargin{2};
%% ====================== End of Checking Inputs ======================= %%

%% ========================== Main Program ============================= %%
fid = fopen(textFile,'wt');
for i = 1:size(cad2print,1)
    fprintf(fid,'%s\n',cad2print(i,:));
    
end
fclose(fid);
%% ======================= End of Main Program ========================= %%

% ---- Outputs
varargout{1} = textFile;
return

function varargout = cell2cad(varargin);
%
% Syntax :
%   cadAll = cell2cad(cellVar, delimChar);
%
% This script creates the cad variable (lines of strings) from a
% bidimensional cell array.
%
% Input Parameters:
%         cellVar          : Matlab cellarray.
%         delimChar        : Delimiter.
%
% Output Parameters:
%         cadAll           : Matrix of characters.
%
%
% Related references:
%
%
% See also:
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% CHUV
% March 14th 2018
% Version $1.0

%% ============================= Checking Inputs ======================= %%


%% ========================== Checking Inputs ========================== %%
cellVar = varargin{1};
if nargin >2
    error('To many inputs');
    return;
end
if nargout >1
    error('To many outputs');
    return;
end

if nargin ~= 2
    delimChar = ',';
else 
    delimChar = varargin{2};
end
%% ====================== End of Checking Inputs ======================= %%

%% ========================== Main Program ============================= %%
cadAll = '';
for i = 1:size(cellVar,1)
    cadLine = '';
    for j = 1:size(cellVar,2)
        tempVar = cellVar{i,j};
        if isnumeric(tempVar);
            tempVar = num2str(tempVar);
        elseif ischar(tempVar)
            
        end
        cadLine = [cadLine delimChar tempVar];
    end
    cadAll = strvcat(cadAll,cadLine);
end
cadAll(:,1) = [];
%% ======================= End of Main Program ========================= %%

% ---- Outputs
varargout{1} = cadAll;
return;



