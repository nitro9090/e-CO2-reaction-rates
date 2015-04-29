function [ output_args ] = Exporttxt(filename, Names, variable, VariableValue,data)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


[~,NameSize] = size(Names);
exportdata(1,2:NameSize+1) = Names;
exportdata(1,1) = cellstr(variable);
exportdata2(:,1) = VariableValue;
exportdata2(:,2:NameSize+1) = data(:,1:NameSize);
[datasize,~] = size(VariableValue);
exportdata(2:datasize+1,:) = num2cell(exportdata2);

%filedata = dataset(exportdata)
dlmcell(filename, exportdata);

end

