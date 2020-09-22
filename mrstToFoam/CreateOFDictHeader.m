function [] = CreateOFDictHeader(fileID,class,location,object)
%% CreateOFfileHeader 
%  This function writes an OF dictionary header with optional additional
%  infos.

fprintf(fileID,'FoamFile\n');
fprintf(fileID,'{\n');
fprintf(fileID,'    version     2.0;\n');
fprintf(fileID,'    format      ascii;\n');
fprintf(fileID,'    class       %s;\n',class);
fprintf(fileID,'    location    "%s";\n',location);
fprintf(fileID,'    object      %s;\n',object);
fprintf(fileID,'}\n');
end

