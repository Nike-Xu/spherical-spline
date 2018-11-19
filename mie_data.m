function [thetagrid,pdata]=mie_data(siz,wavl,realp,imagp,nang_hf)
% generate some spherical data from Mie theory 
% input: 
% 1. siz : particle radius in um (could be array)
% 2. wavl: wavelength in um 
% 3. realp : real part of the refractive index 
% 4. img : imaginary part of the rafractive index 
% 5. nang_hf:     angle number from 0~pi/2;
% output : 
% 1. thetagrid: mesh data of phi and theta  
% 2. pdata: phase function data 

addpath('./mie')  ;
sizex=siz*2.0*pi/wavl;                 % size parameters; 
refin=realp+imagp*1i ;                 % refractive index; 
angs=2*(nang_hf-1)+1;                  % total angles; 
NS = length(sizex);                    % number of sizes ;  
QE(NS)=zeros; QS(NS)=zeros;GF(NS)=zeros; EC(NS)=zeros; 
QB(NS)=zeros; AB(NS)=zeros; p11(NS,angs)=zeros; 


 for k=1:NS      
    [ss1,ss2,QE(k),QS(k),QB(k),GF(k)]=mie(sizex(k),refin,nang_hf); 
    phf=0.50*(abs(ss2).^2+abs(ss1).^2); 
    EC(k)=QE(k)*pi*sizex(k)^2; 
    AB(k)=QS(k)/QE(k);      
    p11(k,1:end)=phf*4.0*pi/EC(k)/AB(k); 
   
 end
 
 
Angs=0:180/(angs-1):180;
thetagrid=pi*Angs/180.0;
pdata=log(p11); 

return
end