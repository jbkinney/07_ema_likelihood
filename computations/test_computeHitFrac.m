clear all
close all

data = loadMukData('abf1');
dataSet = quantizeData(data, 20);
models = loadEndModelsAndEmats('../results/abf1_muk/', 1000);

model = models(1);
x = runmodelondata(model, dataSet, data, 100000, 1, 2)

