function straightLineMechanism_Ka
% Kinematic analysis of slider-crank mechanism in body-coordinates
% Position, velocity, and accleration analysis from t = 0 to t = 20 s
% Coordinate-Partitioning Method
clc

%% Position analysis
q = [1 1.732 0.41 1.0 -1.732 5.873455 4 1.7321 pi/3 6 1.7321 5*pi/3 4 -1.7321 5*pi/3 6 -1.7321 pi/3 1.5 0 0]';

t = 0:0.001:20;N = size(t,2);% N = number of time steps
pcoordsall = zeros(21,N); % position coordinates at each time step
for i = 1:N
    [phi7,~,~] = driver(t(i));
    q(21) = phi7;
    pcoordsall(:,i) = Newton_Raphson(q,21,1e-6,100,@constraints,t(i));
end


%% Velocity Analysis
qdot = zeros(21,1); % initializing the velocities
vcoordsall  = zeros(21,N); % velocities at each time step
% D=20x21; 
for i = 1:N
   [~,phi7dot,~] = driver(t(i));
   qdot(21) = phi7dot;
   [~,D] = constraints(t(i),pcoordsall(:,i));
   Dnew = D(:,1:20);
   rhs = -D(:,21)*qdot(21);
   qv = Dnew\rhs;
   vcoordsall(:,i) = [qv(1:20,1)' qdot(21)]';
end

%% Acceleration analysis
qddot = zeros(21,1);
acoordsall  = zeros(21,N);
for i = 1:N
   [~,~,phi1ddot] = driver(t(i));
   qddot(21) = phi1ddot;
   [~,D] = constraints(t(i),pcoordsall(:,i));
   Dnew = D(:,1:20);
   rhs = gamma(pcoordsall(:,i),vcoordsall(:,i))-D(:,21)*qddot(21);
   qa = Dnew\rhs;
   acoordsall(:,i) = [qa(1:20,1)' qddot(21)]';
end
l1 = 8.7; l2 = 8.7;
l3 = 4; l4 = 4; l5 = 4; l6 = 4;
l7 = 3;

YE = pcoordsall(17,:)+sin(pcoordsall(18,:))*l6/2;
figure
plot(t,YE)
title('Y position of Joint E')
XE = pcoordsall(16,:)+cos(pcoordsall(18,:))*l6/2;
figure
plot(t,XE)
title('X position of Joint E')

VE = vcoordsall(17,:)+cos(pcoordsall(18,:))*l6/2.*vcoordsall(17,:);

figure
plot(t,VE)
title('Y Vel position of Joint E')

AE = acoordsall(17,:)+cos(pcoordsall(18,:))*l6/2.*(vcoordsall(17,:).^2)-sin(pcoordsall(18,:)).*acoordsall(18,:)*l6/2;
figure
plot(t,AE)
title('Y Acc position of Joint E')


save stkadata.mat t pcoordsall vcoordsall acoordsall

end
function output = gamma(q,qdot)

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

phi1dot = qdot(3);phi2dot = qdot(6);phi3dot = qdot(9);
phi4dot = qdot(12);phi5dot = qdot(15);phi6dot = qdot(18);
phi7dot = qdot(21);

% Relative Positions of Joints with respect to the CM of Links


% Relative Positions of Joints with respect to the CM of Links

s_O_0 = [1.5 0]';
s_A_0 = [-1.5 0]';
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



output = [S_O_7*phi7dot^2;
          S_A_1*phi1dot^2;
          S_A_2*phi2dot^2;
          -S_B_1*phi1dot^2+S_B_3*phi3dot^2;
          -S_B_1*phi1dot^2+S_B_4*phi4dot^2;
          -S_C_7*phi7dot^2+S_C_3*phi3dot^2;
          -S_C_7*phi7dot^2+S_C_5*phi5dot^2;
          -S_D_2*phi2dot^2+S_D_5*phi5dot^2;
          -S_D_2*phi2dot^2+S_D_6*phi6dot^2;
          -S_E_4*phi4dot^2+S_E_6*phi6dot^2];

        
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

s_O_0 = [1.5 0]';
s_A_0 = [-1.5 0]';
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



r_O = [-1.5 0]';
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
function [phi7,phi7dot,phi7ddot] = driver(t)
        
     phi7 = sin(t);
     phi7dot = cos(t);
     phi7ddot = -sin(t);
 end

function output = A(phi)
output = [cos(phi) -sin(phi);sin(phi) cos(phi)];
end
