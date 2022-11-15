clear all, close all, clc
%% Read the files in the directory

folderSrc = '*/ecospeclib-all/*.txt';

matrixFilelist = dir(folderSrc);

% only consider the files showing the spectrum
indices = [];
for n = 1:length(matrixFilelist)
    if ~contains(matrixFilelist(n).name,'spectrum')
        indices = [indices n];
    end
end
matrixFilelist(indices) = [];

%% Extract the useful data

database = [];
for file = 1:size(matrixFilelist,1)
      file_name = matrixFilelist(file).name;
      material = readFile(file_name);
      
%       database = [database material];
end


%%

function [material] = readFile(file_name)

    material = struct;
    values = [];
    
    % open the file in mode read
    fileID = fopen(file_name,'r');
    
    % read the content of the file
    A = fscanf(fileID,'%c');
    
    % close the file for modifications
    fclose(fileID);
    
    % split the content based on line breaks
    tmp = splitlines(A);
    
    % remove the empty cells
    tmp = tmp(~cellfun('isempty',tmp));
    
    % add the field name and value to the structure
    for i = 1:length(tmp)
        B = regexp(tmp(i), ':', 'split','once');
        if length(B{1}) >= 2
            fieldname = strrep( strrep(B{1}{1},' ','_') ,'.','');
            material = setfield(material, fieldname ,B{1}{2});
        else
            C = split(B{1},'	');
            values = [values str2double(C)];
        end
        material = setfield(material,'values',values');
    end
end