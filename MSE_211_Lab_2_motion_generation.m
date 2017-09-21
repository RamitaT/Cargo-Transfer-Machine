function MSE_211_Lab_2_motion_generation
close all; clear all;
% This is the main file where the inputs and type of problem are defined

global inputs 

%Define the Precision Points (Position and Orientation)
%Position Vectors Px and Py (INPUT POSITIONS)
% project
inputs.Px(1) = 5;
inputs.Py(1) = 4;
inputs.Px(2) = 3.5;
inputs.Py(2) = 4;
inputs.Px(3) = 2;
inputs.Py(3) = 3;
inputs.Px(4) = 1;
inputs.Py(4) = 0;

%Define Orientation in radians (INPUT ANGLES)
%project
inputs.ang(1) = 0; %In radians
inputs.ang(2) = -5*pi/180;   %In radians
inputs.ang(3) = -5*pi/180;   %In radians
inputs.ang(4) = 0;   %In radians

% Desired angles with respect to the first angle (DO NOT CHANGE)
inputs.alpha(1)=inputs.ang(2)-inputs.ang(1);
inputs.alpha(2)=inputs.ang(3)-inputs.ang(1);
inputs.alpha(3)=inputs.ang(4)-inputs.ang(1);

% Desired positions with respect to the first positions (DO NOT CHANGE)
inputs.deltax(1)=inputs.Px(2)-inputs.Px(1);
inputs.deltay(1)=inputs.Py(2)-inputs.Py(1);

inputs.deltax(2)=inputs.Px(3)-inputs.Px(1);
inputs.deltay(2)=inputs.Py(3)-inputs.Py(1);

inputs.deltax(3)=inputs.Px(4)-inputs.Px(1);
inputs.deltay(3)=inputs.Py(4)-inputs.Py(1);

motion_generation %Call motion generation function (see below)
end

function motion_generation

%Optimization Left Side
%Initial values (INPUT INITIAL ESTIMATES)
r2x = 0.58;
r2y = 3.42;
r3ax = 2.5;
r3ay = 2.5;
beta(1)= 10*pi/180;  %In radians
beta(2)= -250*pi/180;  %In radians
beta(3)= -250*pi/180;  %In radians

% Initial esimates in vector form (DO NOT CHANGE)
x0=[r2x,r2y,r3ax,r3ay,beta]';

%Call Optimization Algorithm (DO NOT CHANGE)
options = optimoptions(@fminunc,'Algorithm','quasi-newton');
[Lx,LFval]= fminunc(@left_side,x0,options); %Call optimization algorithm

%Optimization Right side
%Initial values (INPUT INITIAL ESTIMATES)
r4x = 0.6;
r4y = 4.4;
r3bx = -2.5;
r3by = 0.5;
gamma(1)= 25*pi/180;  %In radians
gamma(2)= 50*pi/180;  %In radians
gamma(3)= -200*pi/180;  %In radians

% Initial esimates in vector form (DO NOT CHANGE)
x0=[r4x,r4y,r3bx,r3by,gamma]';

%Call Optimization Algorithm (DO NOT CHANGE)
[Rx,RFval]= fminunc(@right_side,x0,options);%Call optimization algorithm

%Call animation (DO NOT CHANGE)
animation(Lx,Rx,LFval,RFval)
end

function F=left_side(x)
global inputs

%Redefine your variables (DO NOT CHANGE)
r2x=x(1);
r2y=x(2);
r3ax=x(3);
r3ay=x(4);
beta(1)=x(5);
beta(2)=x(6);
beta(3)=x(7);

%WRITE THE RELATIVE EQUATIONS FOR EACH PRECISION POINT
f1x = r2x + r3ax + inputs.deltax(1) - r3ax*cos(inputs.alpha(1)) + r3ay*sin(inputs.alpha(1)) - r2x*cos(beta(1)) + r2y*sin(beta(1));
f1y = r2y + r3ay + inputs.deltay(1) - r3ay*cos(inputs.alpha(1)) - r3ax*sin(inputs.alpha(1)) - r2y*cos(beta(1)) - r2x*sin(beta(1));
f2x = r2x + r3ax + inputs.deltax(2) - r3ax*cos(inputs.alpha(2)) + r3ay*sin(inputs.alpha(2)) - r2x*cos(beta(2)) + r2y*sin(beta(2));
f2y = r2y + r3ay + inputs.deltay(2) - r3ay*cos(inputs.alpha(2)) - r3ax*sin(inputs.alpha(2)) - r2y*cos(beta(2)) - r2x*sin(beta(2));
f3x = r2x + r3ax + inputs.deltax(3) - r3ax*cos(inputs.alpha(3)) + r3ay*sin(inputs.alpha(3)) - r2x*cos(beta(3)) + r2y*sin(beta(3));
f3y = r2y + r3ay + inputs.deltay(3) - r3ay*cos(inputs.alpha(3)) - r3ax*sin(inputs.alpha(3)) - r2y*cos(beta(3)) - r2x*sin(beta(3));


%WRITE YOUR OBJECTIVE FUNCTION 
F = sqrt(f1x^2 + f1y^2 + f2x^2 + f2y^2 + f3x^2 + f3y^2);
end

function G=right_side(x)
global inputs

%Redefine your variables (DO NOT CHANGE)
r4x=x(1);
r4y=x(2);
r3bx=x(3);
r3by=x(4);
gamma(1)=x(5);
gamma(2)=x(6);
gamma(3)=x(7);

%WRITE THE RELATIVE EQUATIONS FOR EACH PRECISION POINT
g1x = r4x + r3bx + inputs.deltax(1) - r3bx*cos(inputs.alpha(1)) + r3by*sin(inputs.alpha(1)) - r4x*cos(gamma(1)) + r4y*sin(gamma(1));
g1y = r4y + r3by + inputs.deltay(1) - r3by*cos(inputs.alpha(1)) - r3bx*sin(inputs.alpha(1)) - r4y*cos(gamma(1)) - r4x*sin(gamma(1));
g2x = r4x + r3bx + inputs.deltax(2) - r3bx*cos(inputs.alpha(2)) + r3by*sin(inputs.alpha(2)) - r4x*cos(gamma(2)) + r4y*sin(gamma(2));
g2y = r4y + r3by + inputs.deltay(2) - r3by*cos(inputs.alpha(2)) - r3bx*sin(inputs.alpha(2)) - r4y*cos(gamma(2)) - r4x*sin(gamma(2));
g3x = r4x + r3bx + inputs.deltax(3) - r3bx*cos(inputs.alpha(3)) + r3by*sin(inputs.alpha(3)) - r4x*cos(gamma(3)) + r4y*sin(gamma(3));
g3y = r4y + r3by + inputs.deltay(3) - r3by*cos(inputs.alpha(3)) - r3bx*sin(inputs.alpha(3)) - r4y*cos(gamma(3)) - r4x*sin(gamma(3));

%WRITE YOUR OBJECTIVE FUNCTION 
G = sqrt(g1x^2 + g1y^2 + g2x^2 + g2y^2 + g3x^2 + g3y^2);
end

function animation(Lx,Rx,LFval,RFval)
global inputs 
%In this function you might not have to change anything for the lab 
%(assuming that you obtain a similar solution to mine). However, you will
%probably have to make changes for the project, depending on your outputs.
%THE SECTIONS TO BE CHANGED ARE HIGHLIGHTED WITH LINES / CAPITAL LETTERS

% Precision Points as a vector
P1=[inputs.Px(1),inputs.Py(1)];
P2=[inputs.Px(2),inputs.Py(2)];
P3=[inputs.Px(3),inputs.Py(3)];
P4=[inputs.Px(4),inputs.Py(4)];

% Set up the screen size 
set(0,'Units','pixels');   % Set root units to pixels
dim = get(0,'ScreenSize'); % Get screen size
figure('doublebuffer','on','Position',[0,20,dim(3),dim(4)-100],...  %... line continues in the next row
    'Name','3D Object','NumberTitle','off'); %Create a screen size figure
set(gcf,'color', [1 1 1]) %Background Colour

%______________________________________________________________________
%Drawing the box (THIS PART IS ONLY FOR THE LAB)
box_shape = [0,0; 1,0; 1,1; 0,1; 0,0]; %Creating a 2x2 box
box=NaN(5,2,4); %pre-allocating memory
for i=1:4 % rotating and tranlating the box
box(:,:,i)=[cos(inputs.ang(i))*box_shape(:,1)-sin(inputs.ang(i))*box_shape(:,2),sin(inputs.ang(i))*box_shape(:,1)+cos(inputs.ang(i))*box_shape(:,2)] + repmat([inputs.Px(i), inputs.Py(i)],5,1);
line(box(:,1,i),box(:,2,i)); %ploting the box
end
axis([-4 10 -7 7]); axis equal; hold all;  %CHANGE AXIS SCALE TO FIT YOUR PROBLEM
%______________________________________________________________________

% Redefining the variables of the Left side optimization outputs 
R2=[Lx(1),Lx(2)]; R3a=[Lx(3),Lx(4)]; 
beta(1)=Lx(5); beta(2)=Lx(6); beta(3)=Lx(7);

%Finding the length of the links for the left side
r2=norm(R2); r3a=norm(R3a);

% Redefining the variables of the Right Side optimization outputs 
R4=[Rx(1),Rx(2)]; R3b=[Rx(3),Rx(4)];
gamma(1)=Rx(5); gamma(2)=Rx(6); gamma(3)=Rx(7);

%Finding the length of the links for the right side
r4=norm(R4); r3b=norm(R3b);

%Location of Ground Pivots
GP2 = P1-R2-R3a;
GP4 = P1-R4-R3b;

% Finding the length of the frame and couple links
r1=norm(GP2-GP4);
r3 = norm((R3a-P1)-(R3b-P1)); 

% Defining output angles to be between -pi and pi
for i=1:3
if beta(i)>pi; beta(i)=beta(i)-2*pi; end
if beta(i)<-pi; beta(i)=beta(i)+2*pi; end
if gamma(i)>pi; gamma(i)=gamma(i)-2*pi; end
if gamma(i)<pi; gamma(i)=gamma(i)+2*pi; end
end

% Output script of the solution obtained with the optimization
disp(' ')
fprintf('--------------------------------------------------------------------------------------------------------------\n');
disp(['      LFval           RFval            r1            r2              r3             r4           beta1          beta2          beta3         gamma1       gamma2       gamma3'])
fprintf('--------------------------------------------------------------------------------------------------------------\n');
ffa=[LFval    RFval   r1 r2 r3 r4 beta(1)*180/pi beta(2)*180/pi beta(3)*180/pi gamma(1)*180/pi gamma(2)*180/pi gamma(3)*180/pi];
fprintf('%14.8f %14.8f %14.8f %14.8f %14.8f %14.8f %14.8f %14.8f %14.8f %14.8f %14.8f %14.8f\n',ffa');

% Obtaining input and output angles (theta2 and theta4)
theta2_0= atan2(R2(2),R2(1));
theta2_f=theta2_0+beta(3);
theta4_0= atan2(R4(2),R4(1));
theta4_f=theta4_0+gamma(3);

%Define frame angle (always the same)
theta_1=atan2((GP4(2)-GP2(2)),(GP4(1)-GP2(1)));

%Ploting
%Plot frame link
line([GP2(1),GP4(1)],[GP2(2),GP4(2)], 'LineWidth', 3, 'Color', 'black')
plot(GP2(1),GP2(2), 'o', 'MarkerFaceColor', [0 0 1], 'MarkerSize',8, 'MarkerEdgeColor','k')
plot(GP4(1),GP4(2), 'o', 'MarkerFaceColor', [0 0 1], 'MarkerSize',8, 'MarkerEdgeColor','k')

%Plot r2 in all four configurations
plot([GP2(1);GP2(1)+r2*cos(theta2_0)],[GP2(2);GP2(2)+r2*sin(theta2_0)],'b'); 
plot([GP2(1);GP2(1)+r2*cos(theta2_0+beta(1))],[GP2(2);GP2(2)+r2*sin(theta2_0+beta(1))],'b');
plot([GP2(1);GP2(1)+r2*cos(theta2_0+beta(2))],[GP2(2);GP2(2)+r2*sin(theta2_0+beta(2))],'b');
plot([GP2(1);GP2(1)+r2*cos(theta2_0+beta(3))],[GP2(2);GP2(2)+r2*sin(theta2_0+beta(3))],'b');

%Plot r4 in all four configurations
plot([GP4(1);GP4(1)+r4*cos(theta4_0)],[GP4(2);GP4(2)+r4*sin(theta4_0)],'r'); 
plot([GP4(1);GP4(1)+r4*cos(theta4_0+gamma(1))],[GP4(2);GP4(2)+r4*sin(theta4_0+gamma(1))],'r');
plot([GP4(1);GP4(1)+r4*cos(theta4_0+gamma(2))],[GP4(2);GP4(2)+r4*sin(theta4_0+gamma(2))],'r');
plot([GP4(1);GP4(1)+r4*cos(theta4_0+gamma(3))],[GP4(2);GP4(2)+r4*sin(theta4_0+gamma(3))],'r');

%Plot r3 in all four configurations
plot([GP2(1)+r2*cos(theta2_0);GP4(1)+r4*cos(theta4_0);P1(1);GP2(1)+r2*cos(theta2_0)],[GP2(2)+r2*sin(theta2_0);GP4(2)+r4*sin(theta4_0);P1(2);GP2(2)+r2*sin(theta2_0)],'g')
plot([GP2(1)+r2*cos(theta2_0+beta(1));GP4(1)+r4*cos(theta4_0+gamma(1));P2(1);GP2(1)+r2*cos(theta2_0+beta(1))],[GP2(2)+r2*sin(theta2_0+beta(1));GP4(2)+r4*sin(theta4_0+gamma(1));P2(2);GP2(2)+r2*sin(theta2_0+beta(1))],'g')
plot([GP2(1)+r2*cos(theta2_0+beta(2));GP4(1)+r4*cos(theta4_0+gamma(2));P3(1);GP2(1)+r2*cos(theta2_0+beta(2))],[GP2(2)+r2*sin(theta2_0+beta(2));GP4(2)+r4*sin(theta4_0+gamma(2));P3(2);GP2(2)+r2*sin(theta2_0+beta(2))],'g')
plot([GP2(1)+r2*cos(theta2_0+beta(3));GP4(1)+r4*cos(theta4_0+gamma(3));P4(1);GP2(1)+r2*cos(theta2_0+beta(3))],[GP2(2)+r2*sin(theta2_0+beta(3));GP4(2)+r4*sin(theta4_0+gamma(3));P4(2);GP2(2)+r2*sin(theta2_0+beta(3))],'g')

%__________________________________________________________________________
%Definine range of theta2 (MAKE SURE IT WILL ROTATE IN THE RIGHT DIRECTION)
%The mechanism can rotate clockwise or counterclockwise, also the magnitude 
%of the angle is important for example if you move counterclockwise from 
%theta2_0=100 (deg) to theta2_f= 20 (deg), then you have to add 360 deg to 
%theta2_f, as shown below 
theta_2=theta2_0:(2*pi+theta2_f-theta2_0)/100:2*pi+theta2_f;
%__________________________________________________________________________

%__________________________________________________________________________
%If the coupler plate is inverted in your animation, ADD A NEGATIVE SIGN TO
%THE phi angle, i.e., phi=-acos(...) 
%Defining coupler plate
R3=R3a-R3b;
phi=-acos(dot([R3a(1) R3a(2)],[R3(1) R3(2)])/(r3*r3a));
r3a=norm(R3a);
%__________________________________________________________________________

% Pre-allocating memory 
Theta_3=NaN(length(theta_2),1); Theta_4=NaN(length(theta_2),1); 

%Drawing the staircase
steps = [0,0; 5,0; 5,4; 8,4; 8,4];
line(steps(:,1),steps(:,2), 'LineWidth', 2, 'color', 'black'); axis equal
axis([-4 10 -7 7]); axis equal; hold all;  %CHANGE AXIS SCALE TO FIT YOUR PROBLEM


grid minor
%Animating Mechanism
for i=1:length(theta_2)
    %Displacement Analysis
    a = r1*cos(theta_1) - r2*cos(theta_2(i));
    b = r1*sin(theta_1) - r2*sin(theta_2(i));
    c = (r3^2 - a^2 -r4^2 - b^2)/(2*r4);     
    [fx,theta_4] = NewtRaph(atan2(b,a) + atan2(sqrt(a^2+b^2-c^2),c),10E-6,a,b,c);
    theta_3 = get_theta3(theta_4,a,b);
    
    %Plotting Links
    R(2)=line([GP2(1),GP2(1)+r2*cos(theta_2(i))],[GP2(2),GP2(2)+r2*sin(theta_2(i))], 'LineWidth', 3);
    R(3)=line([GP2(1)+r2*cos(theta_2(i)), GP2(1)+r2*cos(theta_2(i))+r3a*cos(theta_3-phi), GP4(1)+r4*cos(theta_4),GP2(1)+r2*cos(theta_2(i))],[GP2(2)+r2*sin(theta_2(i)),GP2(2)+r2*sin(theta_2(i))+r3a*sin(theta_3-phi), GP4(2)+r4*sin(theta_4),GP2(2)+r2*sin(theta_2(i))], 'LineWidth', 2);
    R(4)=line([GP4(1), GP4(1)+r4*cos(theta_4)],[GP4(2),GP4(2)+r4*sin(theta_4)], 'LineWidth', 3);
    
    %Filling translucent Plate 
    Plate=fill([GP2(1)+r2*cos(theta_2(i)), GP2(1)+r2*cos(theta_2(i))+r3a*cos(theta_3-phi), GP4(1)+r4*cos(theta_4),GP2(1)+r2*cos(theta_2(i))],[GP2(2)+r2*sin(theta_2(i)),GP2(2)+r2*sin(theta_2(i))+r3a*sin(theta_3-phi), GP4(2)+r4*sin(theta_4),GP2(2)+r2*sin(theta_2(i))],'b');
    set(Plate,'facealpha',0.5)
    
    %Plotting joints
    C(1)=plot(GP2(1)+r2*cos(theta_2(i)),GP2(2)+r2*sin(theta_2(i)), 'o', 'MarkerFaceColor', [0 0 1], 'MarkerSize',8, 'MarkerEdgeColor','k');
    C(2)=plot(GP4(1)+r4*cos(theta_4),GP4(2)+r4*sin(theta_4), 'o', 'MarkerFaceColor', [0 0 1], 'MarkerSize',8, 'MarkerEdgeColor','k');
    
    %_________________________________________
    % Scaling and pausing animation SCALE YOUR ANIMATION BASED ON YOUR WORKSPACE 
    axis equal; axis([-4 10 -7 7])
    pause(0.1)
    %_________________________________________
    
    % Turn off links and joints
    if i < length(theta_2)
        for ii=2:4; set(R(ii),'visible','off'); end
        set(C,'visible','off')
        set(Plate,'visible','off')
    end
    
    % Collect all Thetas (In case you want to plot them)
    Theta_3(i)=theta_3;
    Theta_4(i)=theta_4;
end

% close all
% 
% theta2 = theta_2;
% theta_3 = Theta_3;
% theta_4 = Theta_4;
% 
% figure (1)
% title('Angular Displacement of link R3 and R4');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Angular Displacement (deg)'); % y-axis label
% hold on;
% plot(radtodeg(theta2),radtodeg(theta_3),radtodeg(theta2),radtodeg(theta_4))
% legend('Angular Displacement (R3)','Angular Displacement (R4)')
% 
% %Find center of mass
% for i = 1:length(theta2)
%     x2(i) = (r2/2)*cos(theta2(i));
%     y2(i) = (r2/2)*sin(theta2(i));
%     x3(i) = r2*cos(theta2(i))+(r3/2)*cos(theta_3(i));
%     y3(i) = r2*sin(theta2(i))+(r3/2)*sin(theta_3(i));
%     x4(i) = (r4/2)*cos(theta_4(i));
%     y4(i) = (r4/2)*sin(theta_4(i));
% end
% 
% figure (2)
% title('Center of Mass of Link R2');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Position'); % y-axis label
% hold on;
% plot(radtodeg(theta2),x2,radtodeg(theta2),y2)
% legend('x2','y2')
% 
% figure (3)
% title('Center of Mass of Link R3');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Position'); % y-axis label
% hold on;
% plot(radtodeg(theta2),x3,radtodeg(theta2),y3)
% legend('x3','y3')
% 
% figure (4)
% title('Center of Mass of Link R4');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Position'); % y-axis label
% hold on;
% plot(radtodeg(theta2),x4,radtodeg(theta2),y4)
% legend('x4','y4')
% 
% %Analytic Solution
% %Angular Velocity and Angular Acceleration
% theta2_v = 1;    %theta2. = 1 rad/s
% 
% for i = 1:length(theta2)
%    A = [-r3*sin(theta_3(i)) r4*sin(theta_4(i));
%        r3*cos(theta_3(i)) -r4*cos(theta_4(i))];
%    B = [r2*sin(theta2(i))*theta2_v;
%        -r2*cos(theta2(i))*theta2_v];
%    C = A\B;
%    theta3_v(i) = C(1);
%    theta4_v(i) = C(2);
%    D = [r2*cos(theta2(i))*(theta2_v)^2+r3*cos(theta_3(i))*(theta3_v(i))^2-r4*cos(theta_4(i))*(theta4_v(i))^2;
%        r2*sin(theta2(i))*(theta2_v)^2+r3*sin(theta_3(i))*(theta3_v(i))^2-r4*sin(theta_4(i))*(theta4_v(i))^2];
%    E = A\D;
%    theta3_a(i) = E(1);
%    theta4_a(i) = E(2);
% end
% 
% figure (5)
% title('Angular Veloctiy of link R3 and R4 (Analytic Solution)');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Angular Velocity (deg/s)'); % y-axis label
% hold on;
% plot(radtodeg(theta2),radtodeg(theta3_v),radtodeg(theta2),radtodeg(theta4_v))
% legend('Angular Velocity (R3)','Angular Velocity (R4)')
% 
% figure (6)
% title('Angular Acceleration of link R3 and R4 (Analytic Solution)');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Angular Acceleration (deg/s)'); % y-axis label
% hold on;
% plot(radtodeg(theta2),radtodeg(theta3_a),radtodeg(theta2),radtodeg(theta4_a))
% legend('Angular Acceleration (R3)','Angular Acceleration (R4)')
% 
% %Linear Velocity and Linear Acceleration
% for i = 1:length(theta2)
%     v2x(i) = -(r2/2)*theta2_v*sin(theta2(i));
%     v2y(i) = (r2/2)*theta2_v*cos(theta2(i));
%     v3x(i) = -r2*theta2_v*sin(theta2(i))-(r3/2)*theta3_v(i)*sin(theta_3(i));
%     v3y(i) = r2*theta2_v*cos(theta2(i))+(r3/2)*theta3_v(i)*cos(theta_3(i));
%     v4x(i) = -(r4/2)*theta4_v(i)*sin(theta_4(i));
%     v4y(i) = (r4/2)*theta4_v(i)*cos(theta_4(i));
%     a2x(i) = -(r2/2)*cos(theta2(i))*(theta2_v)^2;
%     a2y(i) = -(r2/2)*sin(theta2(i))*(theta2_v)^2;
%     a3x(i) = -r2*cos(theta2(i))*(theta2_v)^2-(r3/2)*sin(theta_3(i))*theta3_a(i)-(r3/2)*cos(theta_3(i))*(theta3_v(i))^2;
%     a3y(i) = -r2*sin(theta2(i))*(theta2_v)^2+(r3/2)*cos(theta_3(i))*theta3_a(i)-(r3/2)*sin(theta_3(i))*(theta3_v(i))^2;
%     a4x(i) = -(r4/2)*sin(theta_4(i))*theta4_a(i)-(r4/2)*cos(theta_4(i))*(theta4_v(i))^2;
%     a4y(i) = (r4/2)*cos(theta_4(i))*theta4_a(i)-(r4/2)*sin(theta_4(i))*(theta4_v(i))^2;
% end
% 
% figure (7)
% title('Linear Velocity of the center of mass of link R2 (Analytic Solution)');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Linear Velocity (m/s)'); % y-axis label
% hold on;
% plot(radtodeg(theta2),v2x,radtodeg(theta2),v2y)
% legend('v2x','v2y')
% 
% figure (8)
% title('Linear Velocity of the center of mass of link R3 (Analytic Solution)');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Linear Velocity (m/s)'); % y-axis label
% hold on;
% plot(radtodeg(theta2),v3x,radtodeg(theta2),v3y)
% legend('v3x','v3y')
% 
% figure (9)
% title('Linear Velocity of the center of mass of link R4 (Analytic Solution)');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Linear Velocity (m/s)'); % y-axis label
% hold on;
% plot(radtodeg(theta2),v4x,radtodeg(theta2),v4y)
% legend('v4x','v4y')
% 
% figure (10)
% title('Linear Acceleration of the center of mass of link R2 (Analytic Solution)');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Linear Acceleration (m/s)'); % y-axis label
% hold on;
% plot(radtodeg(theta2),a2x,radtodeg(theta2),a2y)
% legend('a2x','a2y')
% 
% figure (11)
% title('Linear Acceleration of the center of mass of link R3 (Analytic Solution)');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Linear Acceleration (m/s)'); % y-axis label
% hold on;
% plot(radtodeg(theta2),a3x,radtodeg(theta2),a3y)
% legend('a3x','a3y')
% 
% figure (12)
% title('Linear Acceleration of the center of mass of link R4 (Analytic Solution)');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Linear Acceleration (m/s)'); % y-axis label
% hold on;
% plot(radtodeg(theta2),a4x,radtodeg(theta2),a4y)
% legend('a4x','a4y')
% 
% %Numerical Solution
% %Center Finite-Difference (Basic Version)
% h = pi/36;
% for i = 2:length(theta2)-1
%     theta2_num(i-1) = theta2(i);    %Shifting theta2 element
%     %Angular Velocity
%     theta3_v_num(i-1) = (theta_3(i+1)-theta_3(i-1))/(2*h);
%     theta4_v_num(i-1) = (theta_4(i+1)-theta_4(i-1))/(2*h);
%     %Angular Aceleration
%     theta3_a_num(i-1) = (theta_3(i+1)-2*theta_3(i)+theta_3(i-1))/(h^2);
%     theta4_a_num(i-1) = (theta_4(i+1)-2*theta_4(i)+theta_4(i-1))/(h^2);
%     %Linear Velocity
%     v2x_num(i-1) = (x2(i+1)-x2(i-1))/(2*h);
%     v2y_num(i-1) = (y2(i+1)-y2(i-1))/(2*h);
%     v3x_num(i-1) = (x3(i+1)-x3(i-1))/(2*h);
%     v3y_num(i-1) = (y3(i+1)-y3(i-1))/(2*h);
%     v4x_num(i-1) = (x4(i+1)-x4(i-1))/(2*h);
%     v4y_num(i-1) = (y4(i+1)-y4(i-1))/(2*h);
%     %Linear Acceleration
%     a2x_num(i-1) = (x2(i+1)-2*x2(i)+x2(i-1))/(h^2);
%     a2y_num(i-1) = (y2(i+1)-2*y2(i)+y2(i-1))/(h^2);
%     a3x_num(i-1) = (x3(i+1)-2*x3(i)+x3(i-1))/(h^2);
%     a3y_num(i-1) = (y3(i+1)-2*y3(i)+y3(i-1))/(h^2);
%     a4x_num(i-1) = (x4(i+1)-2*x4(i)+x4(i-1))/(h^2);
%     a4y_num(i-1) = (y4(i+1)-2*y4(i)+y4(i-1))/(h^2);
% end
% 
% figure (13)
% title('Angular Veloctiy of link R3 and R4 (Numerical Solution: Basic Version)');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Angular Velocity (deg/s)'); % y-axis label
% hold on;
% plot(radtodeg(theta2_num),radtodeg(theta3_v_num),radtodeg(theta2_num),radtodeg(theta4_v_num))
% legend('Angular Velocity (R3)','Angular Velocity (R4)')
% 
% figure (14)
% title('Angular Acceleration of link R3 and R4 (Numerical Solution: Basic Version)');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Angular Acceleration (deg/s)'); % y-axis label
% hold on;
% plot(radtodeg(theta2_num),radtodeg(theta3_a_num),radtodeg(theta2_num),radtodeg(theta4_a_num))
% legend('Angular Acceleration (R3)','Angular Acceleration (R4)')
% 
% figure (15)
% title('Linear Velocity of the center of mass of link R2 (Numerical Solution: Basic Version)');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Linear Velocity (m/s)'); % y-axis label
% hold on;
% plot(radtodeg(theta2_num),v2x_num,radtodeg(theta2_num),v2y_num)
% legend('v2x','v2y')
% 
% figure (16)
% title('Linear Velocity of the center of mass of link R3 (Numerical Solution: Basic Version)');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Linear Velocity (m/s)'); % y-axis label
% hold on;
% plot(radtodeg(theta2_num),v3x_num,radtodeg(theta2_num),v3y_num)
% legend('v3x','v3y')
% 
% figure (17)
% title('Linear Velocity of the center of mass of link R4 (Numerical Solution: Basic Version)');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Linear Velocity (m/s)'); % y-axis label
% hold on;
% plot(radtodeg(theta2_num),v4x_num,radtodeg(theta2_num),v4y_num)
% legend('v4x','v4y')
% 
% figure (18)
% title('Linear Acceleration of the center of mass of link R2 (Numerical Solution: Basic Version)');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Linear Acceleration (m/s)'); % y-axis label
% hold on;
% plot(radtodeg(theta2_num),a2x_num,radtodeg(theta2_num),a2y_num)
% legend('a2x','a2y')
% 
% figure (19)
% title('Linear Acceleration of the center of mass of link R3 (Numerical Solution: Basic Version)');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Linear Acceleration (m/s)'); % y-axis label
% hold on;
% plot(radtodeg(theta2_num),a3x_num,radtodeg(theta2_num),a3y_num)
% legend('a3x','a3y')
% 
% figure (20)
% title('Linear Acceleration of the center of mass of link R4 (Numerical Solution: Basic Version)');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Linear Acceleration (m/s)'); % y-axis label
% hold on;
% plot(radtodeg(theta2_num),a4x_num,radtodeg(theta2_num),a4y_num)
% legend('a4x','a4y')
% 
% %Center Finite-Difference (Expanded Version)
% %h = pi/36;
% for i = 3:length(theta2)-2
%     theta2_numE(i-2) = theta2(i);    %Shifting theta2 element
%     %Angular Velocity
%     theta3_v_numE(i-2) = (-theta_3(i+2)+8*theta_3(i+1)-8*theta_3(i-1)+theta_3(i-2))/(12*h);
%     theta4_v_numE(i-2) = (-theta_4(i+2)+8*theta_4(i+1)-8*theta_4(i-1)+theta_4(i-2))/(12*h);
%     %Angular Aceleration
%     theta3_a_numE(i-2) = (-theta_3(i+2)+16*theta_3(i+1)-30*theta_3(i)+16*theta_3(i-1)-theta_3(i-2))/(12*h^2);
%     theta4_a_numE(i-2) = (-theta_4(i+2)+16*theta_4(i+1)-30*theta_4(i)+16*theta_4(i-1)-theta_4(i-2))/(12*h^2);
%     %Linear Velocity
%     v2x_numE(i-2) = (-x2(i+2)+8*x2(i+1)-8*x2(i-1)+x2(i-2))/(12*h);
%     v2y_numE(i-2) = (-y2(i+2)+8*y2(i+1)-8*y2(i-1)+y2(i-2))/(12*h);
%     v3x_numE(i-2) = (-x3(i+2)+8*x3(i+1)-8*x3(i-1)+x3(i-2))/(12*h);
%     v3y_numE(i-2) = (-y3(i+2)+8*y3(i+1)-8*y3(i-1)+y3(i-2))/(12*h);
%     v4x_numE(i-2) = (-x4(i+2)+8*x4(i+1)-8*x4(i-1)+x4(i-2))/(12*h);
%     v4y_numE(i-2) = (-y4(i+2)+8*y4(i+1)-8*y4(i-1)+y4(i-2))/(12*h);
%     %Linear Acceleration
%     a2x_numE(i-2) = (-x2(i+2)+16*x2(i+1)-30*x2(i)+16*x2(i-1)-x2(i-2))/(12*h^2);
%     a2y_numE(i-2) = (-y2(i+2)+16*y2(i+1)-30*y2(i)+16*y2(i-1)-y2(i-2))/(12*h^2);
%     a3x_numE(i-2) = (-x3(i+2)+16*x3(i+1)-30*x3(i)+16*x3(i-1)-x3(i-2))/(12*h^2);
%     a3y_numE(i-2) = (-y3(i+2)+16*y3(i+1)-30*y3(i)+16*y3(i-1)-y3(i-2))/(12*h^2);
%     a4x_numE(i-2) = (-x4(i+2)+16*x4(i+1)-30*x4(i)+16*x4(i-1)-x4(i-2))/(12*h^2);
%     a4y_numE(i-2) = (-y4(i+2)+16*y4(i+1)-30*y4(i)+16*y4(i-1)-y4(i-2))/(12*h^2);
% end
% 
% figure (21)
% title('Angular Veloctiy of link R3 and R4 (Numerical Solution: Expanded Version)');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Angular Velocity (deg/s)'); % y-axis label
% hold on;
% plot(radtodeg(theta2_numE),radtodeg(theta3_v_numE),radtodeg(theta2_numE),radtodeg(theta4_v_numE))
% legend('Angular Velocity (R3)','Angular Velocity (R4)')
% 
% figure (22)
% title('Angular Acceleration of link R3 and R4 (Numerical Solution: Expanded Version)');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Angular Acceleration (deg/s)'); % y-axis label
% hold on;
% plot(radtodeg(theta2_numE),radtodeg(theta3_a_numE),radtodeg(theta2_numE),radtodeg(theta4_a_numE))
% legend('Angular Acceleration (R3)','Angular Acceleration (R4)')
% 
% figure (23)
% title('Linear Velocity of the center of mass of link R2 (Numerical Solution: Expanded Version)');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Linear Velocity (m/s)'); % y-axis label
% hold on;
% plot(radtodeg(theta2_numE),v2x_numE,radtodeg(theta2_numE),v2y_numE)
% legend('v2x','v2y')
% 
% figure (24)
% title('Linear Velocity of the center of mass of link R3 (Numerical Solution: Expanded Version)');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Linear Velocity (m/s)'); % y-axis label
% hold on;
% plot(radtodeg(theta2_numE),v3x_numE,radtodeg(theta2_numE),v3y_numE)
% legend('v3x','v3y')
% 
% figure (25)
% title('Linear Velocity of the center of mass of link R4 (Numerical Solution: Expanded Version)');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Linear Velocity (m/s)'); % y-axis label
% hold on;
% plot(radtodeg(theta2_numE),v4x_numE,radtodeg(theta2_numE),v4y_numE)
% legend('v4x','v4y')
% 
% figure (26)
% title('Linear Acceleration of the center of mass of link R2 (Numerical Solution: Expanded Version)');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Linear Acceleration (m/s)'); % y-axis label
% hold on;
% plot(radtodeg(theta2_numE),a2x_numE,radtodeg(theta2_numE),a2y_numE)
% legend('a2x','a2y')
% 
% figure (27)
% title('Linear Acceleration of the center of mass of link R3 (Numerical Solution: Expanded Version)');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Linear Acceleration (m/s)'); % y-axis label
% hold on;
% plot(radtodeg(theta2_numE),a3x_numE,radtodeg(theta2_numE),a3y_numE)
% legend('a3x','a3y')
% 
% figure (28)
% title('Linear Acceleration of the center of mass of link R4 (Numerical Solution: Expanded Version)');
% xlabel('Theta2 (deg)'); % x-axis label
% ylabel('Linear Acceleration (m/s)'); % y-axis label
% hold on;
% plot(radtodeg(theta2_numE),a4x_numE,radtodeg(theta2_numE),a4y_numE)
% legend('a4x','a4y')
% 
 end


