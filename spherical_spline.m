function [abe,Rc_dense ]=spherical_spline(thetas,phi,sph_data,...
    theta_dense,phi_dense)

% some parameters to be set
deg=0;                               % degree of the spherical harmonics
delta=0.00000005;                    % smoothness; if delta=0, no smoothness
beta=1.0;                           

% data column
Yc=reshape(sph_data,[],1);          % y colume

% G matrix 
thetas_g=thetas-pi/2;               %  adjust to -pi/2 to pi/2
Gmatrix=greenmatrix(thetas_g,phi,thetas_g,phi); 
assignin('base','Gmatrix',Gmatrix);

% B matrix 
Np=size(Gmatrix,1);
Bmt=delta*beta*eye([Np Np]);

% A matrix 
Amt=Amatrix(deg,thetas,phi);    % 
assignin('base','Amatrix',Amt);

% add smoothness to G 
Gmatrix=Gmatrix+Bmt;

% c constants and a coefficients;
cb=(Amt*(Gmatrix\(Amt)'));
ab=(Amt*(Gmatrix\Yc));
c_const=cb\ab;
a_cof=Gmatrix\(Amt)'*c_const - (Gmatrix\Yc);

assignin('base','c_const',c_const);
assignin('base','a_cof',a_cof);

% construct the interpolation 
Rc_sph=Amt'*c_const-Gmatrix*a_cof;

% error estimation
abse=sum((Rc_sph-Yc).^2);
abe=sqrt(abse);

% compute G matrix dense 
the_d=theta_dense-pi/2.0;               % -pi/2 to pi/2
Gmatrix_dense=greenmatrix(thetas_g,phi,the_d,phi_dense); 
assignin('base','Gmatrix_dense',Gmatrix_dense);

% compute A matrix dense
Amt_dense=Amatrix(deg,theta_dense,phi_dense);  % 
assignin('base','Amatrix_dense',Amt_dense);

% dense data construction; 
Rc_dense=Amt_dense'*c_const-Gmatrix_dense*a_cof;

return

end