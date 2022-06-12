mA= 7483.60952; % in MeV
ma=1876.1350;
mb= 1876.1350;
mB= 7483.60952;
Ta=66.80;   % Energy of the beam (Heavy particle) MeV
Q=0; 
z=1; % atomic number of incident particle
Z=14; % atomic number of material through which particle is passing . in this I am taking 14 for Silicon
I1=((9.76+58.8*(Z)^-1.19)*Z)*10^-6; % Ionization potential of material in MeV
for theta=10:1:28
 Eb =(((sqrt(mA*ma*Ta))*cosd(theta)+(sqrt(mA*ma*Ta*(cosd(theta))^2+(mB+mb)*(mB*Q+(mB-mA)*Ta))))/(mB+ma))^2;
 Pb=sqrt(2*Eb*mb);
 EB=Ta-Eb;
 PB=sqrt(2*EB*mB);
 phi=asind((Pb/PB)*sind(theta));
 
    B=0.01; % in cm  this is the thickness of material you want to use


    H=B/cosd(theta); % in cm with change in theta , particle will transverse differnt thickness , so will take hypot as thicness for differnt angles
    dx=H-B+0.01;   % in cm no need to change here  except if you want to change 0.01
   gamma=(ma+Eb)/ma;             % ma and Eb needed to be changed for differnt particles, like for light particle use Eb and ma , for heavy particle use mA and EB
   beta=sqrt(1-(1/gamma)^2);
 
 
 y=(2*Eb*0.511/ma); % same as line 16 comment
 eloss=  ((4*3.14*6.022*10^23*1.44*1.44*10^-29*z*z*Z)/(A*y))*(log(2*y/I1));  % Bohr energy loss
 eloss2= ((10^-3*2*0.3006*0.511*Z*z*z)/(beta*beta*A))*(log(2*gamma*gamma*beta*beta*0.511/I1)-beta*beta); % Bathe Bloch energy loss
  dE=dx*eloss2*rho*10^3;  % dE(MeV)= cm*(Mev*cm^2/mg)*10^3*mg/cm^3; % energy loss after passing through material
 Erem1=Eb-dE;  % it will calculate remaining energy after coming out from material
  fid=fopen('data9x.txt','a+');
                fprintf(fid,'%10.1f %10.4f %10.4f %10.4f\n',Eb',dE',phi',EB');
  fclose(fid);
end
M= load('data9x.txt');
x= M(:,1)';
y1=M(:,2)';
y2=M(:,3)';
y3=M(:,4)';
plot(x,y1,'b','Markersize',2,'LineWidth',2)
hold on
plot(y2,y3,'r','Markersize',2,'LineWidth',2)
%hold on
set(gca,'FontSize',16,'FontWeight','bold','linewidth',2)
r=xlabel('{\bf Energy(light particle)MeV }','interpreter','latex');
s=ylabel('{\bf Q value measured using reaction kinematics(MeV)}','interpreter','latex');
set(r,'FontSize',16);
set(r,'FontWeight','bold');
set(s,'FontSize',16);
set(s,'FontWeight','bold');
 
 
 
 