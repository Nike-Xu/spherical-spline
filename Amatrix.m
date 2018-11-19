function Amn = Amatrix (m,  theta_s,phi)

 % m is the order of the fundamental system 
 % N is the number of data 
 % etaN is the direction data storing 
 % has the form of (theta phi), theta is in (0,pi)
 % phi is in (0, 2pi)
 
[L2,T2]=meshgrid(phi,theta_s);
thetac=reshape(T2,[],1);
phic=reshape(L2,[],1);
etaN=[thetac phic];
 N=size(etaN,1);
 
if (m==0)
    Amn(1,1:N)=sqrt(1/(4.0*pi));
    return
end

if (m~=0)
    Amn((m+1)^2,N)=zeros; 
    Amn(1,1:N)=sqrt(1/(4.0*pi));
    for k=1:m   
     Amn(k^2+1:(k+1)^2,1:N)=sphmonics(k,etaN(:,1),etaN(:,2),N);       
    end
end

% rewrite the function using column matrix 
     function sphm=sphmonics(deg,thetaco,phico,nd)
      Plm = legendre(deg,cos(thetaco)); % adjust the range of theta
      Sinf(deg,nd)=zeros; 
      Cosf(deg,nd)=zeros;
      Al=1:deg; 
      Al=Al';
      Sinf(1:deg,1:nd)=(-1).^(Al).*sin(Al.*phico');
      Cosf(1:deg,1:nd)=(-1).^(Al).*cos(Al.*phico');
      lmf=(2*deg+1)*(factorial(deg-Al))/(4.0*pi)./(factorial(deg+Al));
      lmf=sqrt(2)*sqrt(lmf);
      Norf=repmat(lmf,[1 nd]);
      Pmp=Norf.*Cosf.*Plm(2:deg+1,:); % m>0;
      Pmn=Norf.*Sinf.*Plm(2:deg+1,:); % m<0;
      Pmn=flip(Pmn);                  % m<0;
      norf0=sqrt((2*deg+1)/(4.0*pi)); % m=0;
      Pm0=norf0*Plm(1,:);  
      sphm=[Pmn;Pm0;Pmp];  
      return
     end
  %size(Amn)
 return

end





