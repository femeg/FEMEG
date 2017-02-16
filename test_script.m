%% Test example

% In this example we compute the relative error (RE) as a function of the 
% source eccentricity in a multilayered spherical domain using a low 
% resolution mesh and the first order FEM. We do so for dipolar and 
% quadrupolar sources. For details on the construction of the mesh we refer to:
%
%   + Beltrachini, L., "A finite element solution of the EEG forward problem 
%     for multipolar sources", Submitted

%   This file is part of the FEMEG toolbox.
%   Author: Leandro Beltrachini <BeltrachiniL at cardiff.ac.uk>


%% Prepare

clc,clear

% Load anisotropic spherical model
load Spherical_mesh

% Eccentricity values to be evaluated
Nec=16; 
rec=linspace(0.001,r(end)-0.001,Nec);

% Define dipolar moment
q=[1,0,0]*1e-8; % Tangentially-oriented dipolar source [A.m]
Q=[2,0,0,-1,0,-1]*1e-10; % Tangentially-oriented linear quadrupolar source [A.m^2]

% Choose points on the scalp for computing the RE
ind_sc=sqrt(sum(p.^2,2))>r(1)-1e-5;

% Compute stiffness matrix and preconditioners
M_fo = femeg_stiffness( p, t, D);
[L_fo,U_fo]=ilu(M_fo);


%% Iteration over eccentricity values - Dipolar source 

parfor kk=1:Nec
    
    % Source position
    pos = [0,0,rec(kk)];
    
    % Compute numerical solution
    [b_fo,uinf] = femeg_indep_fs(p,t,pos,q,D,4); % Source vector
    [u_n,flag_n] = qmr(M_fo,b_fo,1e-10,4000,L_fo,U_fo); % Solve the system
    u_n=u_n+uinf;
    u_fsc=u_n-mean(u_n(ind_sc));u_fsc=u_fsc(ind_sc); % Average reference
   
    % Compute analytical solution
    ua_dsc=femeg_anisphere(p(ind_sc,:),pos,q,r,cond);
    ua_med=mean(ua_dsc);ua_dsc=ua_dsc-repmat(ua_med,size(ua_dsc,1),1);

    % Compute errors
    re_fo(kk)=norm(u_fsc-ua_dsc)/norm(ua_dsc);

    disp(['Dipolar source (tangential): ', num2str(kk),',   flag (fo): ', num2str(flag_n)])    
    
end

figure,semilogy(rec/r(end),abs(re_fo),'rs-'),title('Dipolar source (tangential)')
ylim([1e-2,2e-1]),ylabel('Relative error'),xlabel('eccentricity'),grid on


%% Iteration over eccentricity values - Quadrupolar source 

parfor kk=1:Nec
    
    % Source position
    pos = [0,0,rec(kk)];
    
    % Compute numerical solution
    [b_fo,uinf] = femeg_indep_fs(p,t,pos,Q,D,8); % Source vector
    [u_n,flag_n] = qmr(M_fo,b_fo,1e-10,4000,L_fo,U_fo); % Solve the system
    u_n=u_n+uinf;
    u_fsc=u_n-mean(u_n(ind_sc));u_fsc=u_fsc(ind_sc); % Average reference
   
    % Compute analytical solution
    ua_dsc=femeg_anisphere(p(ind_sc,:),pos,Q,r,cond);
    ua_med=mean(ua_dsc);ua_dsc=ua_dsc-repmat(ua_med,size(ua_dsc,1),1);

    % Compute errors
    re_fo(kk)=norm(u_fsc-ua_dsc)/norm(ua_dsc);

    disp(['Linear quadrupolar source (tangential): ', num2str(kk),',   flag (fo): ', num2str(flag_n)])    
    
end

figure,semilogy(rec/r(end),abs(re_fo),'rs-'),title('Linear Quadrupolar source (tangential)')
ylim([3e-2,1e0]),ylabel('Relative error'),xlabel('eccentricity'),grid on
