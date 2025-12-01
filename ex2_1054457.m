%% Structural Dynamics for Aerospace Structures
%  Project02
%  Student: Nikolaos Rigatos (1054457)
%  05/02/2021

%% Clear
clear 
close all
clc

%% Beam Geometry
L = 0.4;            %Beam Length (m)
w = 0.03;           %Beam Width  (m)
h = 0.002;          %Beam Height (m)
A = w*h;            %Beam Cross Section Area    (m^2)
I = w*(h^3)/12;     %Beam Cross Section Inertia (m^4)

%% Material Properties
E = 147e9;          %Modulus of Elasticity (Pa)
v = 0.275;          %Poisson Ratio
G = E/(2*(1+v));    %Modulus of Rigidity   (Pa)
dens = 1578;        %Density               (kg/m^3)

%% Elemental Calculations
N_elements = 12;                            %Number of elements
Node_dofs = 2;                              %Degrees of freedom in each node (w,theta)
Le = L/N_elements;                          %Element Length (m)
N_nodes = N_elements + 1;                   %Total Nodes Number
Total_dofs = N_nodes*Node_dofs;             %Total Degrees of Freedom
N_active_nodes = N_nodes - 2;               %Total active nodes
Active_dofs = N_active_nodes*Node_dofs;     %Total degrees of freedom of actived nodes

%% Discretize the Elemental Space
x_value = 0:N_elements;
x_space = Le.*x_value;

%% Consistent Mass Matrix
%Consistent Mass Matrix for Timoshenko Beams according to Professor...
%Saravanos lecture notes(Structural Dynamics - saam)
mbar_a = dens*A*Le;
mbar_d = dens*I*Le;

Me = [mbar_a/3    0           mbar_a/6     0;
      0           mbar_d/3    0            mbar_d/6;
      mbar_a/6    0           mbar_a/3     0;
      0           mbar_d/6    0            mbar_d/3];
     
%% Initialize Total Consistent Mass Matrix Mtotal (i.e set the size of it)
Mtotal = zeros(Total_dofs,Total_dofs);
     
%% Assembly Total Consistent Mass Matrix
start = 1;
stop = 4;
for noe = 1:N_elements
    Mtotal(start:stop,start:stop) = Mtotal(start:stop,start:stop) + Me;
    start = start + Node_dofs;
    stop = stop + Node_dofs;
end

clear start stop

%Verification of Consistent Mass Matrix
mass_check = sum(sum(Mtotal));
M = dens*L*w*h;
error = abs(mass_check-M)/M;

if (error < 0.01)
    disp("The Mass Matrix is correct");  
end 

%% Define the Elemental Stiffness Matrix
%Reduced integration stiffness matrix for 2-node bending(shear)  according to Professor...
%Saravanos lecture notes (FE2-saam)
Ke = [G*A/Le         -G*A/2           -G*A/Le          -G*A/2;
     -G*A/2      (G*A*Le/4+E*I/Le)     G*A/2      (G*A*Le/4-E*I/Le);
     -G*A/Le          G*A/2            G*A/Le           G*A/2;
     -G*A/2      (G*A*Le/4-E*I/Le)     G*A/2      (G*A*Le/4+E*I/Le)];
      
%% Initialize Total Stiffness Matrix Ktotal (i.e set the size of it)
Ktotal = zeros(Total_dofs,Total_dofs);

%% Assembly Total Stiffness Matrix (without considering boundary conditions)
start = 1;
stop = 4;
for noe = 1:N_elements
    Ktotal(start:stop,start:stop) = Ktotal(start:stop,start:stop) + Ke;
    start = start + Node_dofs;
    stop = stop + Node_dofs;
end

clear start stop
mass_check = sum(sum(Ktotal));

%% Define the Elemental Damping Matrix 
a1 = 0.5/100; 
Ce = a1*Ke;        %Proportional Damping

%% Initialize Total Damping Matrix Ctotal (i.e set the size of it)
Ctotal = zeros(Total_dofs,Total_dofs);

%% Assembly Total Damping Matrix (without considering boundary conditions)
start = 1;
stop = 4;
for noe = 1:N_elements
    Ctotal(start:stop,start:stop) = Ctotal(start:stop,start:stop) + Ce;
    start = start + Node_dofs;
    stop = stop + Node_dofs;
end

clear start stop

 %% Boundary_Conditions
 %Removal of n-1 row and column from matrixes: w_(n-1,n-1)=0
    Mtotal(Total_dofs-1,:) = [];
    Mtotal(:,Total_dofs-1) = [];

    Ktotal(Total_dofs-1,:) = [];
    Ktotal(:,Total_dofs-1) = [];
    
    Ctotal(Total_dofs-1,:) = [];
    Ctotal(:,Total_dofs-1) = [];

 %Removal of 1st,2nd row and column from matrixes
    Mtotal(1:2,:) = [];
    Mtotal(:,1:2) = [];

    Ktotal(1:2,:) = [];
    Ktotal(:,1:2) = [];
    
    Ctotal(1:2,:) = [];
    Ctotal(:,1:2) = [];
    
 %% Natural Frequencies of the system (Hz)
 [V,D] = eig(Ktotal,Mtotal);                     %Calculate the eigenvalues and eigenfrequencies
 [D_sorted] = sort(diag(D), 'ascend');           %Sort the eigenfrequencies from lowest to highest
 natfreq = sqrt(D_sorted)/(2*pi);                %fi = 2p/wi
 ftransverse = natfreq(1:2:end-1,:);
 frotation = natfreq(2:2:end,:);
 frotation = [frotation; natfreq(size(natfreq,1),:)];
   
%% First 5 modes
nmodes = 5;

%Normalization of modes 
Vnorm = zeros(size(V));
for i = 1:size(V,2)
    
           Vnorm(:,i) = V(:,i)/max(abs(V(:,i)));
           
end

Vtransverse = V(1:2:end-1,:);
Vrotation = V(2:2:end,:);
Vrotation = [Vrotation; V(size(V,1),:)];

for i = 1:size(Vtransverse,2)

Vtransverse(:,i) = Vtransverse(:,i)/max(abs(Vtransverse(:,i)));
Vrotation(:,i) = Vrotation(:,i)/max(abs(Vrotation(:,i)));

end

Vtransverse = [zeros(1,size(Vtransverse,2));Vtransverse;zeros(1,size(Vtransverse,2))];
Vrotation = [zeros(1,size(Vrotation,2));Vrotation];

%% Superposition
nmodes3 = 3;

omega_avg = (sqrt(D_sorted(1,1))+sqrt(D_sorted(2,1)))/2;
T = 2*pi/omega_avg;

force = 0.1;
F = zeros(size(Mtotal,1),1);
F(floor(size(F,1)/2),1) = force;         %Set the Force at 1/2 of the beam

mbar = V'*Mtotal*V;                      %M matrix in modal space
kbar = V'*Ktotal*V;                      %K matrix in modal space
fbar = V'*F;                             %F in modal space
cbar = a1*kbar;                          %C in modal space

%% Harmonic imput F = 0.1sin(omega_avg*t)
qa = zeros(size(Mtotal,1),1);
for i=1:size(Mtotal,1)
    
    qa(i,1) = fbar(i,1)/(D_sorted(i,1)*((1-omega_avg^2/D_sorted(i,1))));
    
end

Qa = zeros(nmodes3,length(linspace(0,2.5*T,300)));
t = linspace(0,2.5*T,300);

for i = 1:nmodes3
        
        Qa(i,:) = qa(i,1)*sin(omega_avg*t);
        
end

%Using 1 mode

w_1a = Qa(1,:).*V(floor((3/4)*size(V,1)),1);
theta_1a = Qa(1,:).*V(ceil((3/4)*size(V,1)),1);

%Using 3 modes
w_3a = zeros(1,size(t,2));
theta_3a = zeros(1,size(t,2));
    for i = 1:nmodes3
        
        w_3a = w_3a + Qa(i,:).*V(floor((3/4)*size(V,1)),i);
        theta_3a = theta_3a + Qa(i,:).*V(ceil((3/4)*size(V,1)),i);
        
    end

%% Impulse imput F = 0.1*delta
qb = zeros(size(Mtotal,1),1);
for i=1:size(Mtotal,1)
    
    qb(i,1) = fbar(i,1)/(sqrt(D_sorted(i,1)));
    
end

Qb = zeros(nmodes3,length(linspace(0,2.5*T,300)));
t = linspace(0,2.5*T,300);

for i = 1:nmodes3
        
        Qb(i,:) = qb(i,1)*sin(sqrt(D_sorted(i,1))*t);
        
end

%Using 1 mode
w_1b = Qb(1,:).*V(floor((3/4)*size(V,1)),1);
theta_1b = Qb(1,:).*V(ceil((3/4)*size(V,1)),1);

%Using 3 modes
w_3b = zeros(1,size(t,2));
theta_3b = zeros(1,size(t,2));

    for i = 1:nmodes3
        
        w_3b = w_3b + Qb(i,:).*V(floor((3/4)*size(V,1)),i);
        theta_3b = theta_3b + Qb(i,:).*V(ceil((3/4)*size(V,1)),i);
        
    end
    
%% Non conservative linear MDOF solver for natural frequencies

A11 = zeros(size(Mtotal,1),size(Mtotal,2));
A12 = eye(size(Mtotal,1),size(Mtotal,2));
A21 = -Mtotal\Ktotal;
A22 = -Mtotal\Ctotal;
A_ss= [A11, A12; A21, A22];

[Vd,Dd,Wd] = eig(A_ss,'vector');     %Dd are the eigenfrequencies, Vd the right eigenvectors and Wd the left eigenvectors
Dsorted_d = flip(Dd);
natdamped = abs(Dsorted_d(1:2:end))/(2*pi);
ftransv_d = natdamped(1:2:end-1,:);
frot_d = natdamped(2:2:end,:);
frot_d = [frot_d; natdamped(size(natdamped,1),:)];

%% Poles and left/right eigenvectors
Poles = Dsorted_d;
righteigen = sort(Vd,'ascend');
lefteigen = sort(Wd,'ascend');

%% Transient Dynamic Response
T1 = (2*pi)/sqrt(D_sorted(1));
Tmax = (2*pi)/sqrt(D_sorted(end));

dt1 = T1/5;
dt2 = Tmax/5;
tmax = 0.05;

t1 = 0:dt1:tmax;
t2 = 0:dt2:tmax;


%uminus = u(i-1), u = ui and uplus = u(n+1)
uminus  = zeros(size(Mtotal,1),1,length(t1));         
u       = zeros(size(Mtotal,1),1,length(t1));
uplus   = zeros(size(Mtotal,1),1,length(t1));

%Constants for dt1
k1 = Mtotal/dt1^2;
b1 = Ktotal - 2*Mtotal/dt1^2;
c1 = Mtotal/dt1^2;
%Constants for dt2
k2 = Mtotal/dt2^2;
b2 = Ktotal - 2*Mtotal/dt2^2;
c2 = Mtotal/dt2^2;

%Boundary Conditions
u0 = 0;
uminus(:,:,1) = u0;
u(:,:,1) = u0;

P = zeros(size(Mtotal,1),1,length(t1));
t_1=0;
udt1 = zeros(1,length(t1));

for i = 1:length(t1)
       
       P(floor(size(Mtotal,1)),1,i) = force*sin(omega_avg*t_1);         %Calculate the force for each time step
       uplus(:,:,i) = k1\(P(:,:,i) - b1*u(:,:,i) - c1*uminus(:,:,i));   %u(i+1) = (P-b1*ui-c1*u(i-1)/K
       udt1(i) = uplus(floor((3/4)*size(Mtotal,1)),1,i);                %Calculate the response at 3/4 of beam
       
       u(:,:,i+1) = uplus(:,:,i);                                       %For next loop
       uminus(:,:,i+1) = u(:,:,i);                                      %ui = u(i+1), u(i+1)=u(i)
       
       t_1 = t_1 + dt1;

end

uminus  = zeros(size(Mtotal,1),1,length(t1));
u       = zeros(size(Mtotal,1),1,length(t1));
uplus   = zeros(size(Mtotal,1),1,length(t1));

u0 = 0;
uminus(:,:,1) = u0;
u(:,:,1) = u0;


P = zeros(size(Mtotal,1),1,length(t2));
t_2=0;

udt2 = zeros(1,length(t1));

for i = 1:length(t2)
       
       P(floor(size(Mtotal,1)),1,i) = force*sin(omega_avg*t_2);
       uplus(:,:,i) = k2\(P(:,:,i) - b2*u(:,:,i) - c2*uminus(:,:,i));
       udt2(i) = uplus(floor((3/4)*size(Mtotal,1)),1,i);
       
       u(:,:,i+1) = uplus(:,:,i);
       uminus(:,:,i+1) = u(:,:,i);
       
       t_2 = t_2 + dt2;

end

%% Excel Exports
%Question 1: Timoshkeno Consistent Mass Matrix
folder = [pwd, '/Question1'];
    if ~exist(folder, 'dir')
        mkdir(folder);   
    end
  FileName = 'Timoshenko Consistent Mass Matrix.xlsx';
    fullFileName = fullfile(folder,FileName);
    xlswrite(fullFileName,Mtotal);
    
%Question 2: Natural Frequencies of the System
    folder = [pwd, '/Question2'];
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
    FileName = 'Natural Frequencies of System';
    fullFileName = fullfile(folder,FileName);
    xlswrite(fullFileName,natfreq);
    
%Question 5: Damping Matrix
  folder = [pwd, '/Question5'];
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
    FileName = 'Damping Matrix';
    fullFileName = fullfile(folder,FileName);
    xlswrite(fullFileName,Ctotal);

%Question 7: Poles and left/right eigenvectors
  folder = [pwd, '/Question7'];
    if ~exist(folder, 'dir')
        mkdir(folder);
    end 
    FileName = 'Right Eigenvectors';
    fullFileName = fullfile(folder,FileName);
    xlswrite(fullFileName,righteigen);

    FileName = 'Left Eigenvectors';
    fullFileName = fullfile(folder,FileName);
    xlswrite(fullFileName,lefteigen);
 
%% Figure Exports
%Question 2: Natural Frequencies of 1)System 2)Transverse 3)Rotational
fig1 = figure();
bar(natfreq,'k')
title(('System Natural Frequencies'),'FontSize',14);
xlabel('N frequency','FontSize',14)
ylabel('Natural Frequency f (Hz) ','FontSize',14)
saveas(fig1, fullfile([pwd '/Question2'],'System Natural Frequencies.jpg'));
box on

fig2 = figure();
bar(ftransverse,'r')
title(('System Natural Transverse Frequencies'),'FontSize',14);
xlabel('N Transverse Frequency','FontSize',14)
ylabel('Natural Frequency f (Hz) ','FontSize',14)
saveas(fig2, fullfile([pwd '/Question2'],'System Natural Transverse Frequencies.jpg'));
box on

fig3 = figure();
bar(frotation,'b')
title(('System Natural Rotation Frequencies'),'FontSize',14);
xlabel('N Rotation Frequency','FontSize',14)
ylabel('Natural Frequency f (Hz) ','FontSize',14)
saveas(fig3, fullfile([pwd '/Question2'],'System Natural Rotational Frequencies.jpg'));
box on
    
%Question 3: 5 Modes
    folder = [pwd, '/Question3'];
        if ~exist(folder, 'dir')
            mkdir(folder);
        end

        for i = 1:nmodes
    
            fig4 = figure();
            plot(Vtransverse(:,i), '-r')
            grid on
            hold on
            xlabel('Node (transverse)','FontSize',14)
            ylabel('Normalized Transverse Displacement','FontSize',14)
            title(sprintf('Transverse Displacement Mode Shape %d',i),'FontSize',14)
            axis tight
            saveas(fig4, fullfile([pwd '/Question3'],['Transverse Mode Shape' num2str(i) '.jpg']));

            fig5 = figure();
            plot(Vrotation(:,i), '-b')
            grid on
            hold on
            xlabel('Node (rotation)','FontSize',14)
            ylabel('Normalized Rotation','FontSize',14)
            title(sprintf('Rotation Mode Shape %d',i),'FontSize',14)
            axis tight
            saveas(fig5, fullfile([pwd '/Question3'],['Rotational Mode Shape' num2str(i) '.jpg']));

        end
%Question 4

folder = [pwd, '/Question4'];
        if ~exist(folder, 'dir')
            mkdir(folder);
        end
%a

        fig6 = figure();
        hold on
        title('Transverse Deflection of the beam at 3L/4(Harmonic Input)','FontSize',10 );
        xlabel('t [sec]','FontSize',14)
        ylabel('w [m]','FontSize',14)
        xlim([0 2.5*T])
        plot(t,w_1a,'-b')
        plot(t,w_3a,'-r')
        axis tight
        grid on
        legend('1 Mode Shape','3 Mode Shapes','Location','Northeast')
        saveas(fig6, fullfile([pwd '/Question4'],['w_harmonic_superposition.jpg']));
   
        fig7 = figure();
        hold on
        title('Rotation of the beam at 3L/4(Harmonic Input)','FontSize',10 );
        xlabel('t [sec]','FontSize',14)
        ylabel('theta [rad]','FontSize',14)
        xlim([0 2.5*T])
        plot(t,theta_1a,'-b')
        plot(t,theta_3a,'-r')
        axis tight
        grid on
        legend('1 Mode Shape','3 Mode Shapes','Location','Northeast')
        saveas(fig7, fullfile([pwd '/Question4'],['theta_harmonic_superposition.jpg']));
        
        %b
        
        fig8 = figure();
        hold on
        title('Transverse Deflection of the beam at 3L/4(Impulse input)','FontSize',10 );
        xlabel('t [sec]','FontSize',14)
        ylabel('w [m]','FontSize',14)
        xlim([0 2.5*T])
        plot(t,w_1b,'-b')
        plot(t,w_3b,'-r')
        axis tight
        grid on
        legend('1 Mode Shape','3 Mode Shapes','Location','Northeast')
        saveas(fig8, fullfile([pwd '/Question4'],['w_impulse_superposition.jpg']));
        
        fig9 = figure();
        hold on
        title('Rotation of the beam at 3L/4(Impulse Input)','FontSize',10 );
        xlabel('t [sec]','FontSize',14)
        ylabel('theta [rad]','FontSize',14)
        xlim([0 2.5*T])
        plot(t,theta_1b,'-b')
        plot(t,theta_3b,'-r')
        axis tight
        grid on
        legend('1 Mode Shape','3 Mode Shapes','Location','Northeast')
        saveas(fig9, fullfile([pwd '/Question4'],['theta_impulse_superposition.jpg']));

%Question 6
        folder = [pwd, '/Question6'];
        if ~exist(folder, 'dir')
        mkdir(folder);
        end

fig10 = figure();
bar(natdamped,'k')
title(('Damped System Natural Frequencies'),'FontSize',12);
xlabel('N frequency','FontSize',14)
ylabel('Damped Natural Frequency f (Hz) ','FontSize',14)
saveas(fig10, fullfile([pwd '/Question6'],'Damped System Natural Frequencies.jpg'));
box on

fig11 = figure();
bar(ftransv_d,'r')
title(('Damped System Natural Transverse Frequencies'),'FontSize',12);
xlabel('N Transverse Frequency','FontSize',14)
ylabel('Damped Natural Frequency f (Hz) ','FontSize',14)
saveas(fig11, fullfile([pwd '/Question6'],'Damped System Natural Transverse Frequencies.jpg'));
box on

fig12 = figure();
bar(frot_d,'b')
title(('Damped System Natural Rotation Frequencies'),'FontSize',12);
xlabel('N Rotation Frequency','FontSize',14)
ylabel('Damped Natural Frequency f (Hz) ','FontSize',14)
saveas(fig12, fullfile([pwd '/Question6'],'Damped System Natural Rotational Frequencies.jpg'));
box on

fig13 = figure();
hold on
bar(natfreq,'m')
bar(natdamped,'k')
title(('Damped vs Non Damped System Natural Frequencies'),'FontSize',10);
xlabel('N Rotation Frequency','FontSize',14)
ylabel('Natural Frequency f (Hz) ','FontSize',10)
legend('Non Damped Natural Frequency','Damped Natural Frequency','Location','Northwest')
saveas(fig13, fullfile([pwd '/Question6'],'Nat Freq Comparison.jpg'));
box on

%Question 7
        folder = [pwd, '/Question7'];
        if ~exist(folder, 'dir')
        mkdir(folder);
        end
        
fig14 = figure();       
plot(real(Poles), imag(Poles), '*')
plot(Poles, 'o')
title(('Damped system poles'),'FontSize',14);
xlabel('Real','FontSize',14)
ylabel('Im','FontSize',14)
saveas(fig14, fullfile([pwd '/Question7'],'Poles.jpg'));
grid on

%Question 8
 folder = [pwd, '/Question8'];
        if ~exist(folder, 'dir')
        mkdir(folder);
        end
        
fig15 = figure();
        hold on
        plot(t1,udt1,'-r')
        xlabel('t [sec]','FontSize',12)
        ylabel('Transient Response [m]','FontSize',12)
        title(['Transient Response (F = 0.1sin(omega*t) for dt1 = ',num2str(dt1)],'FontSize',9)
        grid on
        saveas(fig15, fullfile([pwd '/Question8'],['Transient Response for dt1.jpg']));
         

fig16 = figure();
        hold on
        plot(t2,udt2,'-r')
        xlabel('t [sec]','FontSize',12)
        ylabel('Transient Response [m]','FontSize',12)
        title(['Transient Response (F = 0.1sin(omega*t) for dt2 = ',num2str(dt2)],'FontSize',9)
        grid on
        saveas(fig16, fullfile([pwd '/Question8'],['Transient Response for dt2.jpg']));
        

