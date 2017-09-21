function theta3 = get_theta3(theta4,a,b)
     r4 = 3.70;
    cos_theta3 = (a + r4*cos(theta4));
    sin_theta3 = (b + r4*sin(theta4));
    theta3 = atan2(sin_theta3,cos_theta3);
end

% theta_3(i) = atan2(b+r4*sin(theta_4(i)),a+r4*cos(theta_4(i)));