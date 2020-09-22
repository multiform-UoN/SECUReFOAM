function [] = mrstGrdeclToPolyMesh(G)
%% Converts mrst processed grdecl grids to polyMesh format

%% Init lists and compute auxiliary stuff

display("Computing auxiliary quantities");
vertNum = G.nodes.num;
vertList = G.nodes.coords;

faceNum = G.faces.num;
faceNodeNum = diff(G.faces.nodePos);
faceNodePos = G.faces.nodePos;
faceNodeId = G.faces.nodes;

owner =G.faces.neighbors(:,1);
neighbor =G.faces.neighbors(:,2);


nofBFaces = numel(G.faces.neighbors) - nnz(G.faces.neighbors);
totFaces = length(owner);

%% Fix owner/neighbour and face direction to conform to OpenFOAM
% first fix owner (it can not have zeroes)
% It is important that boundary faces go last
% A mapping vector that allows reshuffling is created 

display("Converting data");

boundId = 1;
inId = 1;
swapDir = zeros(1,totFaces);
bound = zeros(1,nofBFaces);
internal = zeros(1,totFaces-nofBFaces);

for i=1:length(owner)
    
    % Put in ownerBound if boundary
    if owner(i) == 0 % Boundary, so no neighbour
        owner(i) = neighbor(i);
        bound(boundId) = i;   
        swapDir(i) = 1;
        boundId = 1+boundId;
    elseif neighbor(i) == 0 %Boundary, so no neighbour
        bound(boundId) = i;   
        boundId = 1+boundId;
    elseif owner(i)> neighbor(i)
        internal(inId) = i;
        swapDir(i) = 1;
        inId = inId+1;
    elseif owner(i)< neighbor(i)
        swapDir(i) = 0;
        internal(inId) = i;
        inId = inId+1;
    end       
    
end

%Merge lists
faceMap = [internal bound ];

%% Create polyMesh folder

mkdir polyMesh;

%% Create points file

display("writing points");

fileID = fopen('polyMesh/points','w');
CreateOFDictHeader(fileID,"vectorField","constant/polyMesh","points");
fprintf(fileID,'%i\n(',vertNum);

for i=1:vertNum
    fprintf(fileID,'   (%f %f %f)\n',vertList(i,1),vertList(i,2),vertList(i,3));
end
fprintf(fileID,');');
fclose(fileID);

%% Create faces file

display("Writing faces");

fileID = fopen('polyMesh/faces','w');
CreateOFDictHeader(fileID,"faceList","constant/polyMesh","faces");
fprintf(fileID,'%i\n(\n',faceNum);

%Write faces first
for k=1:totFaces
    i = faceMap(k); %This is the id of the  face with boundary faces last
    
    faceId = faceNodePos(i); %This is the id in the list of face nodes
    fprintf(fileID,'   %i(',faceNodeNum(i));
    
    % Loop over all the nodes in face i
    for j=1:faceNodeNum(i)
        
        if swapDir(i) == 0
            id = faceId+j-1;
        else
            id = faceId + faceNodeNum(i) -j ;
        end
        
        fprintf(fileID,'%i ',faceNodeId(id)-1);
    end
    fprintf(fileID,')\n');
end

fprintf(fileID,');');

fclose(fileID);


%% Create owner file

display("Writing owner");

fileID = fopen('polyMesh/owner','w');
CreateOFDictHeader(fileID,"labelList","constant/polyMesh","owner");
fprintf(fileID,'%i\n(\n',length(owner));

for k=1:length(owner)
    i = faceMap(k); %Again, boundary last
    fprintf(fileID,'    %i \n',owner(i)-1);
end
fprintf(fileID,');');

fclose(fileID);

%% create neighbour file

display("Writing neighbour");

fileID = fopen('polyMesh/neighbour','w');
CreateOFDictHeader(fileID,"labelList","constant/polyMesh","neighbour");
fprintf(fileID,'%i\n(\n',totFaces-nofBFaces);

%Only internal
for k=1:(totFaces-nofBFaces)
    i = internal(k); %Map to internal faces
    fprintf(fileID,'    %i \n',neighbor(i)-1);
end
fprintf(fileID,');');

fclose(fileID);

%% Finally create boundary file

display("Writing boundary");

fileID = fopen('polyMesh/boundary','w');
CreateOFDictHeader(fileID,"polyBoundaryMesh","constant/polyMesh","boundary");

% Only default faces are available now
fprintf(fileID,'1\n(\n');
fprintf(fileID,'    defaultFaces\n');
fprintf(fileID,'    {\n');
fprintf(fileID,'        type            patch;\n');
fprintf(fileID,'        inGroups        List<word> 1(patch);\n');
fprintf(fileID,'        nFaces          %i;\n',nofBFaces);
fprintf(fileID,'        startFace       %i;\n',totFaces-nofBFaces);
fprintf(fileID,'    }\n');
fprintf(fileID,')\n');

fclose(fileID);
end

