% quick code to map the Outcrop Density averages 

load('OutcropDensitiesFull.mat'); 
load('../BlockyLayerBudget1to438.mat', 'isINS'); 

[x y full] = size(Gbottom);
wintMedGb = zeros(x,y,6);
wintModM = zeros(x,y,6);
wb= 35;
we = 63;
for ye=1:6
	myG = Gbottom(:,:,((ye-1)*73+wb):((ye-1)*73+we));
	wintMedGb(:,:,ye) = median(myG,3);
	Gmod = myG;
	Gmod(Gmod<24)==0;
	Gmod(Gmod>=24)==1;
       	wintModM(:,:,ye) = sum(myG.*Gmod,3)./sum(Gmod,3);
		
end
save('WintModMGbot.mat', 'wintMGb','wb' ,'we') 
