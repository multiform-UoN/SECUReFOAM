close all;
clear all;
%% Read/Create grid

display("Reading grdecl file");

%grdeclFile = fullfile(ROOTDIR,  'borzecin.GRDECL');
%grdeclFile = fullfile(ROOTDIR,  'Sleipner_Reference_Model.GRDECL');

mrstModule add deckformat libgeometry 

% Read in grid
sl_file = fullfile(ROOTDIR,  'borzecinNoFlag.GRDECL');
grdecl = readGRDECL(sl_file);

% Second loading of Sleipner Eclispe grid, to get MAPAXES
%fn      = fopen(sl_file);
%gr  = readGRID(fn, fileparts(sl_file), initializeDeck());
% this grdecl contains: GRID, and others. grdecl.GRID contains
% MAPUNITS, MAPAXES, cartDims, COORD, ZCORN, ACTNUM
%fclose(fn);


% Add data loaded from first loading of Sleipner Eclispe grid
%grdecl.MAPAXES = gr.GRID.MAPAXES;
%grdecl.MAPUNITS = gr.GRID.MAPUNITS;
%clear gr sl_file

% Transform coordinates

% If required, recompute X and Y coordinates in terms of
% the provided axes (depths, Z, do not require any
% recomputation)
%coords        = reshape(grdecl.COORD,3,[])';
%coords(:,1:2) = mapAxes(coords(:,1:2), grdecl.MAPAXES);
%coords        = coords';
%grdecl.COORD  = coords(:); clear coords


%% Next, we process the grid and compute geometry

G = processGRDECL(grdecl); 
%G = mcomputeGeometry(G);

%rock        = grdecl2Rock(grdecl, G.cells.indexMap);
%rock.perm   = convertFrom(rock.perm, milli*darcy);

%grdecl         = readGRDECL(grdeclFile);
%grdecl = simpleGrdecl([8, 5, 7],0.1)
%G = processGRDECL(grdecl);

% If true only fields are processed
fieldsOnly = true;

if fieldsOnly == false
    mrstGrdeclToPolyMesh(G)
end

%% Now it is time to write the poro and perm fields

%Create 0 folder
mkdir 0;

%% Write porosity

display("Writing porosity field");

fileID = fopen('0/poro.orig','w');
CreateOFDictHeader(fileID,"volScalarField","0","poro");

fprintf(fileID,'dimensions      [0 0 0 0 0 0 0];\n');
fprintf(fileID,'internalField   nonuniform List<scalar>\n');
fprintf(fileID,'%i\n(\n',length(G.cells.indexMap));

for i=1:length(G.cells.indexMap)
    id = G.cells.indexMap(i);
    fprintf(fileID,'    %f\n',grdecl.PORO(id));
end
fprintf(fileID,');\n');

fprintf(fileID,'boundaryField\n');
fprintf(fileID,'{\n');
fprintf(fileID,'    defaultFaces\n');
fprintf(fileID,'    {\n');
fprintf(fileID,'        type            zeroGradient;\n');
fprintf(fileID,'    }\n');
fprintf(fileID,'}\n');



%% Write permeability

display("Writing permeability field");

fileID = fopen('0/K.orig','w');
CreateOFDictHeader(fileID,"volTensorField","0","K");

fprintf(fileID,'dimensions      [0 2 0 0 0 0 0];\n');
fprintf(fileID,'internalField   nonuniform List<tensor>\n');
fprintf(fileID,'%i\n(\n',length(G.cells.indexMap));

for i=1:length(G.cells.indexMap)
    id = G.cells.indexMap(i);
    
    % Settings for Borzecin
    permxy = grdecl.PERMXY(id)*9.869233e-16;
    permz = permxy*0.1;
    fprintf(fileID,'    ( %e 0 0 0 %e 0 0 0 %e )\n',permxy,permxy,permz);
end
fprintf(fileID,');\n');

fprintf(fileID,'boundaryField\n');
fprintf(fileID,'{\n');
fprintf(fileID,'    defaultFaces\n');
fprintf(fileID,'    {\n');
fprintf(fileID,'        type            zeroGradient;\n');
fprintf(fileID,'    }\n');
fprintf(fileID,'}\n');


