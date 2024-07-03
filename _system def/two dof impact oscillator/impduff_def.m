%% Problem definition of 2 DoF impact Duffing oscillator
% Governing equation of motion
%  m1*x1'' + b1*x1' + b2*(x1'-x2') + k1*x1 + kn1*x1^3 + k2*(x1-x2) 
%    + kn2*(x1-x2)^3 = F*cos(Omega*t)
%  m2*x2'' - b2*(x1'-x2') - k2*(x1-x2) - kn2*(x1-x2)^3 = 0
% Dimensionless form
%  x1'' + 2*zeta*x1' + 2*xi*mu*(x1'-x2') + x1 + eta1*x1^3 + 
%    + fi^2*mu*(x1-x2) + eta2*mu*(x1-x2)^3 = cos(omega*t)
%  x2'' - 2*xi*(x1'-x2') - fi^2*(x1-x2) - eta2*(x1-x2)^3 = 0
% system parameters: 
%   1) omega = Omega/om_n1:     dimensionless excitation frequency
%   2) delta1 = cl1/x_u:        dimensionless clearance (upwards)
%   3) delta2 = cl2/x_u:        dimensionless clearance (downwards)
%   4) mu = m2/m1:              mass ratio
%   5) r:                       coefficient of restitution
%   6) zeta = c1/(2*m1*om_n1):  standalone damping ration
%   7) xi = c2/(2*m2*om_n1):    TMD damping ratio
%   8) fi = om_n2/om_n1:        TMD frequency ratio
%   9) eta1 = (kn1*xu^2)/(m1*om_n1^2):   standalone nonlinear cubic stiffnes term
%   10) eta2 = (kn2*xu^2)/(m2*om_n1^2):   TMD nonlinear cubic stiffnes term
%   with x_u = F/k1, tu = 1/om_n1, om_n1 = sqrt(m1/k1), om_n2 = sqrt(m2/k2);
% Impact events
% condtions: 
%   delta1  = x2 - x1 (impact at upper surface)
%   delta2  = x1 - x2 (impact at lower surface)
% event map:
%    vs = (x1'+mu*x2')/(1+mu)
%    v1 = (r+1)*vs - r*x1'
%    v2 = (r+1)*vs - r*x2'
% Autonomous extension with cos(omega*t), sin(omega*t)
% x1' = x3
% x2' = x4
% x3' = -x1 - fi^2*mu*(x1-x2) - 2*zeta*x(3) - 2*xi*mu*(x3-x4) 
%   - eta1*x1^3 - eta2*mu*(x1-x2)^3 + x5
% x4' = fi^2*(x1-x2) + 2*xi*(x3-x4) + eta2*(x1-x2)^3
% x5' = -omega*x6 + (1-x5^2-x6^2)*x5
% x6' = omega*x5 + (1-x5^2-x6^2)*x6

% Governing vector field
struct.f{1} = @(x,p) [x(3); x(4); -x(1) - p(8)^2*p(4)*(x(1)-x(2)) - 2*p(6)*x(3)...
    - 2*p(7)*p(4)*(x(3)-x(4)) - p(9)*x(1)^3 - p(10)*p(4)*(x(1)-x(2))^3 + x(5);...
    p(8)^2*(x(1)-x(2)) + 2*p(7)*(x(3)-x(4)) + p(10)*(x(1)-x(2))^3;...
    -p(1)*x(6) + (1-x(5)^2-x(6)^2)*x(5); ...
    p(1)*x(5) + (1-x(5)^2-x(6)^2)*x(6)];

% Vector field Jacobian
struct.Ju{1} = @(x,p) [0 0 1 0 0 0; 0 0 0 1 0 0;...
    -(1+p(8)^2*p(4))-3*p(9)*x(1)^2-3*p(10)*p(4)*(x(1)-x(2))^2 ...
    p(8)^2*p(4)+3*p(10)*p(4)*(x(1)-x(2))^2 -(2*p(6)+2*p(7)*p(4)) ...
    2*p(7)*p(4) 1 0;...
    p(8)^2+3*p(10)*(x(1)-x(2))^2 -p(8)^2-3*p(10)*(x(1)-x(2))^2 ...
    2*p(7) -2*p(7) 0 0;...
    0 0 0 0 1-3*x(5)^2-x(6)^2 -p(1)-2*x(5)*x(6);...
    0 0 0 0 p(1)-2*x(5)*x(6) 1-x(5)^2-3*x(6)^2];

% Parameter Jacobian of the vector field
struct.Jp{1} = @(x,p) [zeros(2,10); zeros(2,3) ...
    [-p(8)^2*(x(1)-x(2))-2*p(7)*(x(3)-x(4))-p(10)*(x(1)-x(2))^3 ...
    0 -2*x(3) 2*p(4)*(x(3)-x(4))  -2*p(8)*p(4)*(x(1)-x(2)) -x(1)^3 ...
    -p(4)*(x(1)-x(2))^3;...
    0 0 0 2*(x(3)-x(4)) 2*p(8)*p(4)*(x(1)-x(2)) 0 (x(1)-x(2))^3];
    [-x(6); x(5)] zeros(2,9)];

% Event surfaces (2 normal and 2 grazing events)
struct.h{1} = @(x,p) x(2)-x(1)-p(2);
struct.h{2} = @(x,p) x(2)-x(1)+p(3);
struct.h{3} = struct.h{1};
struct.h{4} = struct.h{2};

% Event maps
struct.g{1} = @(x,p)  [x(1); x(2);...
    (p(5)+1)/(1+p(4))*(x(3)+p(4)*x(4))-p(5)*x(3);...
    (p(5)+1)/(1+p(4))*(x(3)+p(4)*x(4))-p(5)*x(4);
    x(5); x(6)];
struct.g{2} = struct.g{1};
struct.g{3} = @(x,p) x;
struct.g{4} = struct.g{3};

% Mode changes
struct.pm{1} = [1 1];
struct.pm{2} = struct.pm{1};
struct.pm{3} = struct.pm{1};
struct.pm{4} = struct.pm{1};

% Event surface Jacobians
struct.dh{1} = @(x,p) [-1 1 0 0 0 0];
struct.dh{2} = struct.dh{1};
struct.dh{3} = struct.dh{1};
struct.dh{4} = struct.dh{2};

% Event map Jacobians
struct.dg{1} = @(x,p) [eye(2) zeros(2,4);...
    0 0 (p(5)+1)/(1+p(4))-p(5) p(4)*(p(5)+1)/(1+p(4)) 0 0;...
    0 0 (p(5)+1)/(1+p(4)) p(4)*(p(5)+1)/(1+p(4))-p(5) 0 0;...
    zeros(2,4) eye(2)]; 
struct.dg{2} = struct.dg{1};
struct.dg{3} = @(x,p) eye(6);
struct.dg{4} = struct.dg{3};

% Event surface parameter Jacobians
struct.dhp{1} = @(x,p) [0 -1 zeros(1,8)];
struct.dhp{2} = @(x,p) [0  0 1 zeros(1,7)];
struct.dhp{3} = struct.dhp{1};
struct.dhp{4} = struct.dhp{2};

% Event map parameter Jacobian
struct.dgp{1} = @(x,p) [zeros(2,10); zeros(2,3)...
    [-(p(5)+1)/(1+p(4))^2*x(3)+(p(5)+1)/(1+p(4))^2*x(4) ...
    (1/(1+p(4))-1)*x(3)+p(4)/(1+p(4))*x(4); ...
    -(p(5)+1)/(1+p(4))^2*x(3)+(p(5)+1)/(1+p(4))^2*x(4)...
    1/(1+p(4))*x(3)+(p(4)/(1+p(4))-1)*x(4)] zeros(2,5);...
    zeros(2,10)];
struct.dgp{2} = struct.dgp{1};
struct.dgp{3} = @(x,p) zeros(6,10);
struct.dgp{4} = struct.dgp{3};

% system dimensions
struct.dims = [20, 6];
struct.type = 2; % non-smooth ODE

% Plot functions
f_pl{1} = @(t,x) x(1,:);
f_pl{2} = @(t,x) x(2,:)- x(1,:);

