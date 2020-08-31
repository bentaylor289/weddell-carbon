% quick code to analyze the Outcrop Densities... 

load('OutcropDensitiesFull.mat'); 
load('../BlockyLayerBudget1to438.mat', 'isINS'); 

myEdges = [0.1 20 25:0.2:27.4 27.5:0.1:28 28.05:0.05:28.2 28.24 28.27 28.3:0.05: 28.7]; 

numbin = length(myEdges)-1; 
histOut = zeros(438,numbin);
GbottomN = Gbottom(:,78:398,:);
for k=1:438
     h = histogram(GbottomN(:,:,k).*isINS, myEdges);
     histOut(k,:) = h.Values;
    end
save('HistOutcrop1.mat', 'myEdges', 'histOut') 
