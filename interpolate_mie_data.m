
% test data from Mie thoery
%% Clear 
clear ;close all;clc; 
%% generate Mie data;

p_size=20.0 ;  % size in um ; 
wavelg=0.50 ; % in um; 
realp= 1.20;  % real part of the refractive index ; 
imagp = 0.003 ; % imaginary part of the refractive index ; 
numb_haf=600; % number of angles for half sphere; 

[thetas, phfs]=mie_data(p_size, wavelg, realp, imagp,numb_haf);
 phi0=0.5-pi;
 phi=0.5;

theta_flip=-flip(thetas(2:length(thetas)));
%phfs_flip=flip(phfs(2:length(phfs)));

%new_theta=[theta_flip thetas];
new_theta=thetas; 
%new_phi=[phi0 phi];
new_phi=phi0;
%new_phfs=[phfs_flip phfs];
new_phfs=phfs;

% data preperation
nphi=length(new_phi);
sphere_data=repmat(new_phfs,[nphi 1]);
sphere_data=sphere_data';

% dense distribution of phis 
theta_dense=0:(pi/6000):pi;
%%  spline approximations

%theta_dense=thetas;
[abe, RC]=spherical_spline(new_theta, new_phi, sphere_data,...
    theta_dense,new_phi);
Rc_sph=RC';


%% plot some figures 


angles=rad2deg(new_theta);
angles_dense=rad2deg(theta_dense);
figure(1);
plot(log(angles_dense), Rc_sph,'k');
hold on;
plot(log(angles), sphere_data,'r');
legend('fiting','original');
Rc_cubic=spline(new_theta,sphere_data, theta_dense);
figure(2);
plot(log(angles_dense), Rc_cubic,'k');
hold on;
plot(log(angles), sphere_data,'r');
legend('fiting','original');
figure(3);
plot(log(angles_dense), Rc_cubic,'k');
hold on; 
plot(log(angles_dense), Rc_sph,'r');
legend('cubic spline','spherical spline');
dif=sqrt(sum((Rc_cubic-Rc_sph).^2));
figure(4);
plot(angles, sphere_data,'r');
legend('spherical data');












