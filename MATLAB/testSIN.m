clear all;
close all;
clc;

%% Construct training data (can be multi-dimensional)
data.numS = 200;
xi = sort(rand(data.numS,1)*2*pi);
GT = sin(xi).*(xi<pi) + sin(5*xi).*(xi>=pi);
TRAINING_X = xi;
TRAINING_Y = GT + 0.1*randn(data.numS,1);
TESTING_X = linspace(0,2*pi,200)';

%%% Multi-scale implementation = no need to define a specific scale
IN.gammaList = 10.^[-1,1]; %%% here, two regularization weigths tested

%%% Single scale implementation
% IN.gammaList = 10.^[1];
% IN.singleScale = 0.25;

IN.xi = TRAINING_X; %%% size = numTrainingSamples x numDims
IN.yi = TRAINING_Y; %%% size = numTrainingSamples x numDims
IN.x = TESTING_X; %%% size = numTestingSamples x numDims

OUT = interpolateMulti(IN);

%% Plotting
figure('units','normalized','position',[0 0 1 1]);
hold on;
plot(IN.xi,IN.yi,'k.');
axis square; grid on;

for i=1:length(IN.gammaList)
    alphaT = (i-1 + eps)/(length(IN.gammaList)-1 + eps);
    plot(IN.x,OUT.OUT{i},'Color',[1-alphaT,0,alphaT]);
end
