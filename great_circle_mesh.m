function [dsig]=great_circle_mesh(theta0,phi0,theta1,phi1)

% apply to theta range -pi/2 to pi; 
% phi range -pi to pi; 
% compute the great distance angle;

[tm0, tm1]=meshgrid(theta0, theta1);
[pm0, pm1]=meshgrid(phi0, phi1);

cth0_c=cos(theta0); 
cth1_c=cos(theta1); 
sth0_c=sin(theta0);
sth1_c=sin(theta1);

cc12_n=cth1_c'*cth0_c; 
ss12_n=sth1_c'*sth0_c; 
cs12_n=cth1_c'*sth0_c; 
sc12_n=sth1_c'*cth0_c; 

% phi matrix ; 
pmm=abs(pm0-pm1);            
cphi=cos(pmm);   
sphi=sin(pmm); 

% outer product ; 

cth1 = cos(tm1); 
cc12_cphi = kron(cphi,cc12_n);
ss12_m = repmat(ss12_n,size(cphi));
sc12_cphi = kron(cphi,sc12_n);
cs12_m = repmat(cs12_n, size(cphi));
cth1_sphi = kron(sphi,cth1);
uper=sqrt(cth1_sphi.^2. + (cs12_m-sc12_cphi).^2.) ;
lower=ss12_m+cc12_cphi; 
sig=atan2(uper,lower);

% cosine of angle 
dsig=cos(sig);

return
end