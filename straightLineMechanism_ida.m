function straightLineMechanism_ida

l1 = 8.7; l2 = 8.7;
l3 = 4; l4 = 4; l5 = 4; l6 = 4;
l7 = 3;
g = 9.81;
m1 = .3; J1 = m1*l1^2/12;
m2 = .3; J2 = m2*l2^2/12;
m3 = .2; J3 = m3*l3^2/12;
m4 = .2; J4 = m4*l4^2/12;
m5 = .2; J5 = m5*l5^2/12;
m6 = .2; J6 = m6*l6^2/12;
m7 = .2; J7 = m7*l7^2/12;

%% Mass Matrix

M = diag([m1 m1 J1 m2 m2 J2 m3 m3 J3 m4 m4 J4 m5 m5 J5 m6 m6 J6 m7 m7 J7]);
h = [0 -m1*g 0 0 -m2*g 0 0 -m3*g 0 0 -m4*g 0 0 -m5*g 0 0 -m6*g 0 0 -m7*g 0]';
load stkadata
N = size(t,2);
torque = zeros(21,N);

for i=1:N
    [Phi,D] = constraints(t(i),pcoordsall(:,i));
    
    D_driver = zeros(1,21);
    D_driver(1,21) = 1;
    D_all = [D' D_driver'];
    RHS = M*acoordsall(:,i)-h;
    lambda = D_all\RHS;
    h_driver = D_driver'*lambda(21);
    torque(:,i) = h_driver;
end
plot(t,torque(21,:))
xlabel('time')
ylabel('torque applied on the crank')
title('Input Torque vs time')

end

function [Phi,D] = constraints(t,q)
l1 = 8.7; l2 = 8.7;
l3 = 4; l4 = 4; l5 = 4; l6 = 4;
l7 = 3;
x1 = q(1);  y1 = q(2);  phi1 = q(3);
x2 = q(4);  y2 = q(5);  phi2 = q(6);
x3 = q(7);  y3 = q(8);  phi3 = q(9);
x4 = q(10); y4 = q(11); phi4 = q(12);
x5 = q(13); y5 = q(14); phi5 = q(15);
x6 = q(16); y6 = q(17); phi6 = q(18);
x7 = q(19); y7 = q(20); phi7 = q(21);


% Relative Positions of Joints with respect to the CM of Links

s_O_0 = [0 0]';
s_A_0 = [-3 0]';
s_A_1 = [-l1/2 0]';
s_A_2 = [-l2/2 0]';
s_B_1 = [l1/2 0]';
s_B_3 = [l3/2 0]';
s_B_4 = [-l4/2 0]';
s_C_3 = [-l3/2 0]';
s_C_5 = [-l5/2 0]';
s_C_7 = [l7/2 0]';
s_D_2 = [l2/2 0]';
s_D_5 = [l5/2 0]';
s_D_6 = [-l6/2 0]';
s_E_4 = [l4/2 0]';
s_E_6 = [l6/2 0]';
s_O_7 = [-l7/2 0]';


% Position of Joints in terms of Global Coordinates 

S_O_0 = A(0)*s_O_0;
S_A_0 = A(0)*s_A_0;
S_A_1 = A(phi1)*s_A_1;
S_A_2 = A(phi2)*s_A_2;
S_B_1 = A(phi1)*s_B_1;
S_B_3 = A(phi3)*s_B_3;
S_B_4 = A(phi4)*s_B_4;
S_C_3 = A(phi3)*s_C_3;
S_C_5 = A(phi5)*s_C_5;
S_C_7 = A(phi7)*s_C_7;
S_D_2 = A(phi2)*s_D_2;
S_D_5 = A(phi5)*s_D_5;
S_D_6 = A(phi6)*s_D_6;
S_E_4 = A(phi4)*s_E_4;
S_E_6 = A(phi6)*s_E_6;
S_O_7 = A(phi7)*s_O_7;


S_A_1r = A(pi/2)*S_A_1;
S_A_2r = A(pi/2)*S_A_2;
S_B_1r = A(pi/2)*S_B_1;
S_B_3r = A(pi/2)*S_B_3;
S_B_4r = A(pi/2)*S_B_4;
S_C_3r = A(pi/2)*S_C_3;
S_C_5r = A(pi/2)*S_C_5;
S_C_7r = A(pi/2)*S_C_7;
S_D_2r = A(pi/2)*S_D_2;
S_D_5r = A(pi/2)*S_D_5;
S_D_6r = A(pi/2)*S_D_6;
S_E_4r = A(pi/2)*S_E_4;
S_E_6r = A(pi/2)*S_E_6;
S_O_7r = A(pi/2)*S_O_7;



r_O = [0 0]';
r_1 = [x1 y1]';
r_2 = [x2 y2]';
r_3 = [x3 y3]';
r_4 = [x4 y4]';
r_5 = [x5 y5]';
r_6 = [x6 y6]';
r_7 = [x7 y7]';



r_O_0 = r_O + S_O_0;
r_A_0 = r_O + S_A_0;
r_A_1 = r_1 + S_A_1;
r_A_2 = r_2 + S_A_2;
r_B_1 = r_1 + S_B_1;
r_B_3 = r_3 + S_B_3;
r_B_4 = r_4 + S_B_4;
r_C_3 = r_3 + S_C_3;
r_C_5 = r_5 + S_C_5;
r_C_7 = r_7 + S_C_7;
r_D_2 = r_2 + S_D_2;
r_D_5 = r_5 + S_D_5;
r_D_6 = r_6 + S_D_6;
r_E_4 = r_4 + S_E_4;
r_E_6 = r_6 + S_E_6;
r_O_7 = r_7 + S_O_7;

Phi = [r_O_7 - r_O_0;
    r_A_1-r_A_0;
    r_A_2-r_A_0;
    r_B_3-r_B_1;
    r_B_4-r_B_1;
    r_C_3-r_C_7;
    r_C_5-r_C_7;
    r_D_5-r_D_2;
    r_D_6-r_D_2;
    r_E_6-r_E_4];

D = [zeros(2,18) eye(2) S_O_7r;    
    eye(2) S_A_1r zeros(2,18);
    zeros(2,3) eye(2) S_A_2r zeros(2,15);
    -eye(2) -S_B_1r zeros(2,3) eye(2) S_B_3r zeros(2,12);
    -eye(2) -S_B_1r zeros(2,6)  eye(2) S_B_4r zeros(2,9);
    zeros(2,6) eye(2) S_C_3r zeros(2,9) -eye(2) -S_C_7r;
    zeros(2,6) zeros(2,6) eye(2) S_C_5r zeros(2,3) -eye(2) -S_C_7r;
    zeros(2,3) -eye(2) -S_D_2r zeros(2,6) eye(2) S_D_5r zeros(2,6);
    zeros(2,3) -eye(2) -S_D_2r zeros(2,9) eye(2) S_D_6r zeros(2,3);
    zeros(2,9) -eye(2) -S_E_4r zeros(2,3) eye(2) S_E_6r zeros(2,3)];
end

function output = A(phi)
output = [cos(phi) -sin(phi);sin(phi) cos(phi)];
end


