function G3=greenmatrix(theta0,phi0,theta1,phi1)

% apply to theta range -pi/2 to pi; 
% phi range -pi to pi; 

dircos=great_circle_mesh(theta0,phi0,theta1,phi1);

% compute the greens functions 
clog=dircos==1.0;
dlog=dircos==-1.0;

clogm=double(clog);  % cos=1.0
dlogm=double(dlog);  % cos=-1.0

Amins=1-dircos;
Aplus=1+dircos;

Nk=20;  % number of term for summation 
xdir=Amins/2.0;
k=1; 
L2=(xdir).^k./((k).^2);
for k=2:Nk
    L2=L2+(xdir).^k./((k).^2);
end

G1=1-log(Amins).*(log(Aplus)-log(2))-L2-log(2)^2+log(2).*log(Aplus);
G2=G1./(4*pi);
G2(isnan(G2))=0;          % replace the NANs with zeros; 
G3=G2+(1.0/(4.0*pi))*clogm+(1.0/(4.0*pi)-pi/24.0)*dlogm; 

return

end