function data = ReadData(filename)

dataMat = readtable(filename);

varNames = {'LoadingCase','F0','F1','F2','F3D','LVDT'};

data = dataMat(:,1:6);
data.Properties.VariableNames = varNames;

end
