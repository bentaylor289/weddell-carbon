%% code to find the flux thru an arbitrary contour. 

%v = ones(499,378);% just needs to be 378 by 498 and non divergent ish.
%u = ones(499,378); 
% this is pretty nondivergent apparently!
%v= repmat(sin(0.01*(1:499)),[378 1])';
%u = repmat(cos(0.01*(1:378)),[499 1]);
v= repmat(sin(0.01*(1:378)),[499 1]); % divergent and makes a big trans! cool! 
u = repmat(cos(0.01*(1:378)),[499 1]);
isINS = zeros(498,377); 
isINS(DICconc>22600)=1;
C = M(:,1:1420); % an arbitrary contour in matlab contour output form. 
le = C(2,1);
xx = squeeze(C(2, 2:le+1));
yy = squeeze(C(1, 2:le+1));
trans = 0;
shoelace = 0;
for k=1:le-1
    xk = floor(min(xx(k+1), xx(k)));
    yk = floor(min(yy(k+1), yy(k)));
    if (xx(k+1)>= xx(k)) % equals case will go to 0, so no worries.
        vk = v(xk, yk);
    else
        vk = v(xk, yk+1); % or tmp = min(yy(k+1), yy(k)); c= frac(c); c*v(xk, yk) + (1-c)*v(xc, yk+1)
    end
    
    if (yy(k+1)>= yy(k))
        uk = u(xk, yk);
    else
        uk = u(xk+1, yk);
    end
    transk = -(xx(k+1) - xx(k))*vk + (yy(k+1)-yy(k))*uk; 
    trans = trans+ transk;
    
    shoek = xx(k)*yy(k+1) - xx(k+1)*yy(k);
    shoelace = shoelace + shoek;
    
    x1 = xx(k)- xk; x2 = xx(k+1) - xk; 
    y1 = yy(k)- yk; y2 = yy(k+1) - yk;
    A = abs((x2-x1)*(y2-y1))/2;
    if (x1<=x2 & y1<=y2)
        Ak = A+1+(y1-1)*x2;
    elseif (x1<=x2 & y1 > y2)
        Ak = A + x1 +y2 - x1*y2;
    elseif (x1 > x2 & y1 <= y2)
        Ak = A+1 - x1*y2;
    else
        Ak = A + 1 +(x2-1)*y1;
    end
    isINS(xk,yk) = Ak;     
end
area= abs(shoelace/2)

