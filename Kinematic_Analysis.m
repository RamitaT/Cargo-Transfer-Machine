clear all
close all
clc

%Declare variables (theta2 from 0 to 2pi with step of pi/36)
r1 = 3.92;
r2 = 3.23;
r3 = 4.09;
r4 = 3.70;
theta2 =  0:pi/36:2*pi;

%Calculate theta3 and theta4 (using Newton-Raphson method)
%the mechanism doesn't jump when x0 = 1 & x0 = 8;
theta_1 = 0.8980;
for i = 1:length(theta2)
    a = r1*cos(theta_1) - r2*cos(theta2(i));
    b = r1*sin(theta_1) - r2*sin(theta2(i));
    c = (r3^2 - a^2 -r4^2 - b^2)/(2*r4);     
    [fx(i),theta_4(i)] = NewtRaph(atan2(b,a) + atan2(sqrt(a^2+b^2-c^2),c),10E-6,a,b,c);
    theta_3(i) = get_theta3(theta_4(i),a,b); 
end

figure (1)
title('Angular Displacement of link R3 and R4');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Angular Displacement (deg)'); % y-axis label
hold on;
plot(radtodeg(theta2),radtodeg(theta_3),radtodeg(theta2),radtodeg(theta_4))
legend('Angular Displacement (R3)','Angular Displacement (R4)')

% %Animation
% x = 0:pi/36:2*pi;
% for i = 1:length(x)
%     x1 = x(i);
%     fx = theta_3(i);
%     fxx = theta_4(i);
%     plot(x1,fx,'o','MarkerFaceColor','b','Markersize',8)
%     plot(x1,fxx,'o','MarkerFaceColor','c','Markersize',8)
%     axis([0 2*pi -10 10])
%     grid on;
%     M(i) = getframe;
% end

%Find center of mass
x2 = (r2/2)*cos(theta2);
y2 = (r2/2)*sin(theta2);
x3 = r2*cos(theta2)+(r3/2)*cos(theta_3);
y3 = r2*sin(theta2)+(r3/2)*sin(theta_3);
x4 = (r4/2)*cos(theta_4);
y4 = (r4/2)*sin(theta_4);

figure (2)
title('Center of Mass of Link R2');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Position'); % y-axis label
hold on;
plot(radtodeg(theta2),x2,radtodeg(theta2),y2)
legend('x2','y2')

figure (3)
title('Center of Mass of Link R3');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Position'); % y-axis label
hold on;
plot(radtodeg(theta2),x3,radtodeg(theta2),y3)
legend('x3','y3')

figure (4)
title('Center of Mass of Link R4');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Position'); % y-axis label
hold on;
plot(radtodeg(theta2),x4,radtodeg(theta2),y4)
legend('x4','y4')

%Analytic Solution
%Angular Velocity and Angular Acceleration
theta2_v = 1;    %theta2. = 1 rad/s

for i = 1:length(theta2)
   A = [-r3*sin(theta_3(i)) r4*sin(theta_4(i));
       r3*cos(theta_3(i)) -r4*cos(theta_4(i))];
   B = [r2*sin(theta2(i))*theta2_v;
       -r2*cos(theta2(i))*theta2_v];
   C = A\B;
   theta3_v(i) = C(1);
   theta4_v(i) = C(2);
   D = [r2*cos(theta2(i))*(theta2_v)^2+r3*cos(theta_3(i))*(theta3_v(i))^2-r4*cos(theta_4(i))*(theta4_v(i))^2;
       r2*sin(theta2(i))*(theta2_v)^2+r3*sin(theta_3(i))*(theta3_v(i))^2-r4*sin(theta_4(i))*(theta4_v(i))^2];
   E = A\D;
   theta3_a(i) = E(1);
   theta4_a(i) = E(2);
end

figure (5)
title('Angular Veloctiy of link R3 and R4 (Analytic Solution)');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Angular Velocity (deg/s)'); % y-axis label
hold on;
plot(radtodeg(theta2),radtodeg(theta3_v),radtodeg(theta2),radtodeg(theta4_v))
legend('Angular Velocity (R3)','Angular Velocity (R4)')

figure (6)
title('Angular Acceleration of link R3 and R4 (Analytic Solution)');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Angular Acceleration (deg/s)'); % y-axis label
hold on;
plot(radtodeg(theta2),radtodeg(theta3_a),radtodeg(theta2),radtodeg(theta4_a))
legend('Angular Acceleration (R3)','Angular Acceleration (R4)')

%Linear Velocity and Linear Acceleration
for i = 1:length(theta2)
    v2x(i) = -(r2/2)*theta2_v*sin(theta2(i));
    v2y(i) = (r2/2)*theta2_v*cos(theta2(i));
    v3x(i) = -r2*theta2_v*sin(theta2(i))-(r3/2)*theta3_v(i)*sin(theta_3(i));
    v3y(i) = r2*theta2_v*cos(theta2(i))+(r3/2)*theta3_v(i)*cos(theta_3(i));
    v4x(i) = -(r4/2)*theta4_v(i)*sin(theta_4(i));
    v4y(i) = (r4/2)*theta4_v(i)*cos(theta_4(i));
    a2x(i) = -(r2/2)*cos(theta2(i))*(theta2_v)^2;
    a2y(i) = -(r2/2)*sin(theta2(i))*(theta2_v)^2;
    a3x(i) = -r2*cos(theta2(i))*(theta2_v)^2-(r3/2)*sin(theta_3(i))*theta3_a(i)-(r3/2)*cos(theta_3(i))*(theta3_v(i))^2;
    a3y(i) = -r2*sin(theta2(i))*(theta2_v)^2+(r3/2)*cos(theta_3(i))*theta3_a(i)-(r3/2)*sin(theta_3(i))*(theta3_v(i))^2;
    a4x(i) = -(r4/2)*sin(theta_4(i))*theta4_a(i)-(r4/2)*cos(theta_4(i))*(theta4_v(i))^2;
    a4y(i) = (r4/2)*cos(theta_4(i))*theta4_a(i)-(r4/2)*sin(theta_4(i))*(theta4_v(i))^2;
end

figure (7)
title('Linear Velocity of the center of mass of link R2 (Analytic Solution)');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Linear Velocity (m/s)'); % y-axis label
hold on;
plot(radtodeg(theta2),v2x,radtodeg(theta2),v2y)
legend('v2x','v2y')

figure (8)
title('Linear Velocity of the center of mass of link R3 (Analytic Solution)');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Linear Velocity (m/s)'); % y-axis label
hold on;
plot(radtodeg(theta2),v3x,radtodeg(theta2),v3y)
legend('v3x','v3y')

figure (9)
title('Linear Velocity of the center of mass of link R4 (Analytic Solution)');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Linear Velocity (m/s)'); % y-axis label
hold on;
plot(radtodeg(theta2),v4x,radtodeg(theta2),v4y)
legend('v4x','v4y')

figure (10)
title('Linear Acceleration of the center of mass of link R2 (Analytic Solution)');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Linear Acceleration (m/s)'); % y-axis label
hold on;
plot(radtodeg(theta2),a2x,radtodeg(theta2),a2y)
legend('a2x','a2y')

figure (11)
title('Linear Acceleration of the center of mass of link R3 (Analytic Solution)');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Linear Acceleration (m/s)'); % y-axis label
hold on;
plot(radtodeg(theta2),a3x,radtodeg(theta2),a3y)
legend('a3x','a3y')

figure (12)
title('Linear Acceleration of the center of mass of link R4 (Analytic Solution)');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Linear Acceleration (m/s)'); % y-axis label
hold on;
plot(radtodeg(theta2),a4x,radtodeg(theta2),a4y)
legend('a4x','a4y')

%Numerical Solution
%Center Finite-Difference (Basic Version)
h = pi/36;
for i = 2:length(theta2)-1
    theta2_num(i-1) = theta2(i);    %Shifting theta2 element
    %Angular Velocity
    theta3_v_num(i-1) = (theta_3(i+1)-theta_3(i-1))/(2*h);
    theta4_v_num(i-1) = (theta_4(i+1)-theta_4(i-1))/(2*h);
    %Angular Aceleration
    theta3_a_num(i-1) = (theta_3(i+1)-2*theta_3(i)+theta_3(i-1))/(h^2);
    theta4_a_num(i-1) = (theta_4(i+1)-2*theta_4(i)+theta_4(i-1))/(h^2);
    %Linear Velocity
    v2x_num(i-1) = (x2(i+1)-x2(i-1))/(2*h);
    v2y_num(i-1) = (y2(i+1)-y2(i-1))/(2*h);
    v3x_num(i-1) = (x3(i+1)-x3(i-1))/(2*h);
    v3y_num(i-1) = (y3(i+1)-y3(i-1))/(2*h);
    v4x_num(i-1) = (x4(i+1)-x4(i-1))/(2*h);
    v4y_num(i-1) = (y4(i+1)-y4(i-1))/(2*h);
    %Linear Acceleration
    a2x_num(i-1) = (x2(i+1)-2*x2(i)+x2(i-1))/(h^2);
    a2y_num(i-1) = (y2(i+1)-2*y2(i)+y2(i-1))/(h^2);
    a3x_num(i-1) = (x3(i+1)-2*x3(i)+x3(i-1))/(h^2);
    a3y_num(i-1) = (y3(i+1)-2*y3(i)+y3(i-1))/(h^2);
    a4x_num(i-1) = (x4(i+1)-2*x4(i)+x4(i-1))/(h^2);
    a4y_num(i-1) = (y4(i+1)-2*y4(i)+y4(i-1))/(h^2);
end

figure (13)
title('Angular Veloctiy of link R3 and R4 (Numerical Solution: Basic Version)');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Angular Velocity (deg/s)'); % y-axis label
hold on;
plot(radtodeg(theta2_num),radtodeg(theta3_v_num),radtodeg(theta2_num),radtodeg(theta4_v_num))
legend('Angular Velocity (R3)','Angular Velocity (R4)')

figure (14)
title('Angular Acceleration of link R3 and R4 (Numerical Solution: Basic Version)');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Angular Acceleration (deg/s)'); % y-axis label
hold on;
plot(radtodeg(theta2_num),radtodeg(theta3_a_num),radtodeg(theta2_num),radtodeg(theta4_a_num))
legend('Angular Acceleration (R3)','Angular Acceleration (R4)')

figure (15)
title('Linear Velocity of the center of mass of link R2 (Numerical Solution: Basic Version)');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Linear Velocity (m/s)'); % y-axis label
hold on;
plot(radtodeg(theta2_num),v2x_num,radtodeg(theta2_num),v2y_num)
legend('v2x','v2y')

figure (16)
title('Linear Velocity of the center of mass of link R3 (Numerical Solution: Basic Version)');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Linear Velocity (m/s)'); % y-axis label
hold on;
plot(radtodeg(theta2_num),v3x_num,radtodeg(theta2_num),v3y_num)
legend('v3x','v3y')

figure (17)
title('Linear Velocity of the center of mass of link R4 (Numerical Solution: Basic Version)');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Linear Velocity (m/s)'); % y-axis label
hold on;
plot(radtodeg(theta2_num),v4x_num,radtodeg(theta2_num),v4y_num)
legend('v4x','v4y')

figure (18)
title('Linear Acceleration of the center of mass of link R2 (Numerical Solution: Basic Version)');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Linear Acceleration (m/s)'); % y-axis label
hold on;
plot(radtodeg(theta2_num),a2x_num,radtodeg(theta2_num),a2y_num)
legend('a2x','a2y')

figure (19)
title('Linear Acceleration of the center of mass of link R3 (Numerical Solution: Basic Version)');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Linear Acceleration (m/s)'); % y-axis label
hold on;
plot(radtodeg(theta2_num),a3x_num,radtodeg(theta2_num),a3y_num)
legend('a3x','a3y')

figure (20)
title('Linear Acceleration of the center of mass of link R4 (Numerical Solution: Basic Version)');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Linear Acceleration (m/s)'); % y-axis label
hold on;
plot(radtodeg(theta2_num),a4x_num,radtodeg(theta2_num),a4y_num)
legend('a4x','a4y')

%Center Finite-Difference (Expanded Version)
%h = pi/36;
for i = 3:length(theta2)-2
    theta2_numE(i-2) = theta2(i);    %Shifting theta2 element
    %Angular Velocity
    theta3_v_numE(i-2) = (-theta_3(i+2)+8*theta_3(i+1)-8*theta_3(i-1)+theta_3(i-2))/(12*h);
    theta4_v_numE(i-2) = (-theta_4(i+2)+8*theta_4(i+1)-8*theta_4(i-1)+theta_4(i-2))/(12*h);
    %Angular Aceleration
    theta3_a_numE(i-2) = (-theta_3(i+2)+16*theta_3(i+1)-30*theta_3(i)+16*theta_3(i-1)-theta_3(i-2))/(12*h^2);
    theta4_a_numE(i-2) = (-theta_4(i+2)+16*theta_4(i+1)-30*theta_4(i)+16*theta_4(i-1)-theta_4(i-2))/(12*h^2);
    %Linear Velocity
    v2x_numE(i-2) = (-x2(i+2)+8*x2(i+1)-8*x2(i-1)+x2(i-2))/(12*h);
    v2y_numE(i-2) = (-y2(i+2)+8*y2(i+1)-8*y2(i-1)+y2(i-2))/(12*h);
    v3x_numE(i-2) = (-x3(i+2)+8*x3(i+1)-8*x3(i-1)+x3(i-2))/(12*h);
    v3y_numE(i-2) = (-y3(i+2)+8*y3(i+1)-8*y3(i-1)+y3(i-2))/(12*h);
    v4x_numE(i-2) = (-x4(i+2)+8*x4(i+1)-8*x4(i-1)+x4(i-2))/(12*h);
    v4y_numE(i-2) = (-y4(i+2)+8*y4(i+1)-8*y4(i-1)+y4(i-2))/(12*h);
    %Linear Acceleration
    a2x_numE(i-2) = (-x2(i+2)+16*x2(i+1)-30*x2(i)+16*x2(i-1)-x2(i-2))/(12*h^2);
    a2y_numE(i-2) = (-y2(i+2)+16*y2(i+1)-30*y2(i)+16*y2(i-1)-y2(i-2))/(12*h^2);
    a3x_numE(i-2) = (-x3(i+2)+16*x3(i+1)-30*x3(i)+16*x3(i-1)-x3(i-2))/(12*h^2);
    a3y_numE(i-2) = (-y3(i+2)+16*y3(i+1)-30*y3(i)+16*y3(i-1)-y3(i-2))/(12*h^2);
    a4x_numE(i-2) = (-x4(i+2)+16*x4(i+1)-30*x4(i)+16*x4(i-1)-x4(i-2))/(12*h^2);
    a4y_numE(i-2) = (-y4(i+2)+16*y4(i+1)-30*y4(i)+16*y4(i-1)-y4(i-2))/(12*h^2);
end

figure (21)
title('Angular Veloctiy of link R3 and R4 (Numerical Solution: Expanded Version)');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Angular Velocity (deg/s)'); % y-axis label
hold on;
plot(radtodeg(theta2_numE),radtodeg(theta3_v_numE),radtodeg(theta2_numE),radtodeg(theta4_v_numE))
legend('Angular Velocity (R3)','Angular Velocity (R4)')

figure (22)
title('Angular Acceleration of link R3 and R4 (Numerical Solution: Expanded Version)');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Angular Acceleration (deg/s)'); % y-axis label
hold on;
plot(radtodeg(theta2_numE),radtodeg(theta3_a_numE),radtodeg(theta2_numE),radtodeg(theta4_a_numE))
legend('Angular Acceleration (R3)','Angular Acceleration (R4)')

figure (23)
title('Linear Velocity of the center of mass of link R2 (Numerical Solution: Expanded Version)');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Linear Velocity (m/s)'); % y-axis label
hold on;
plot(radtodeg(theta2_numE),v2x_numE,radtodeg(theta2_numE),v2y_numE)
legend('v2x','v2y')

figure (24)
title('Linear Velocity of the center of mass of link R3 (Numerical Solution: Expanded Version)');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Linear Velocity (m/s)'); % y-axis label
hold on;
plot(radtodeg(theta2_numE),v3x_numE,radtodeg(theta2_numE),v3y_numE)
legend('v3x','v3y')

figure (25)
title('Linear Velocity of the center of mass of link R4 (Numerical Solution: Expanded Version)');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Linear Velocity (m/s)'); % y-axis label
hold on;
plot(radtodeg(theta2_numE),v4x_numE,radtodeg(theta2_numE),v4y_numE)
legend('v4x','v4y')

figure (26)
title('Linear Acceleration of the center of mass of link R2 (Numerical Solution: Expanded Version)');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Linear Acceleration (m/s)'); % y-axis label
hold on;
plot(radtodeg(theta2_numE),a2x_numE,radtodeg(theta2_numE),a2y_numE)
legend('a2x','a2y')

figure (27)
title('Linear Acceleration of the center of mass of link R3 (Numerical Solution: Expanded Version)');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Linear Acceleration (m/s)'); % y-axis label
hold on;
plot(radtodeg(theta2_numE),a3x_numE,radtodeg(theta2_numE),a3y_numE)
legend('a3x','a3y')

figure (28)
title('Linear Acceleration of the center of mass of link R4 (Numerical Solution: Expanded Version)');
xlabel('Theta2 (deg)'); % x-axis label
ylabel('Linear Acceleration (m/s)'); % y-axis label
hold on;
plot(radtodeg(theta2_numE),a4x_numE,radtodeg(theta2_numE),a4y_numE)
legend('a4x','a4y')