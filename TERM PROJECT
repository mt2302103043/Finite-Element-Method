E = 200e9;I = 30*10^-6;q0 = 20000;P = 50000;R = 10;A = 25;

x = [0,1];
k1 = elestiff(E,I,x);

x = [1,2];
k2 = elestiff(E,I,x);

x = [2,3];
k3 = elestiff(E,I,x);

x = [3,4];
k4 = elestiff(E,I,x);

x = [4,5];
k5 = elestiff(E,I,x);

x = [5,6];
k6 = elestiff(E,I,x);

x = [6,7];
k7 = elestiff(E,I,x);

x = [7,8];
k8 = elestiff(E,I,x);

x = [8,9];
k9 = elestiff(E,I,x);

% ASSEMBlY into Global Matrix

K = zeros(20,20);
F = zeros(20,1);
K(1:4,1:4) = k1(1:4,1:4);
K(3:6,3:6) = K(3:6,3:6) + k2(1:4,1:4);
K(5:8,5:8) = K(5:8,5:8) + k3(1:4,1:4);
K(7:10,7:10) = K(7:10,7:10) + k4(1:4,1:4);
K(9:12,9:12) = K(9:12,9:12) + k5(1:4,1:4);
K(11:14,11:14) = K(11:14,11:14) + k6(1:4,1:4);
K(13:16,13:16) = K(13:16,13:16) + k7(1:4,1:4);
K(15:18,15:18) = K(15:18,15:18) + k8(1:4,1:4);
K(17:20,17:20) = K(17:20,17:20) + k9(1:4,1:4);
F(19) = F(19) + P;


Kreduce = K(3:20,3:20);
Freduce = F(3:20);


ureduce = (Kreduce)\Freduce;

un = [0;0;ureduce];
Freac = K*un;
xn = [0,1,2,3,4,5,6,7,8,9];
xnume = []; unume = [];
for e1 = 1:9
    x_n = xn(e1:e1+1);
    u_n = un((e1-1)*2+1:2*(e1+1));
    le = x_n(2)-x_n(1);
    xi = (-1:0.2:1)';
    Nx = [(1-xi)/2,(1+xi)/2];
    N1 = (2-3*xi+xi.^3)/4;
    N2 = (1-xi-xi.^2+xi.^3)/4;
    N3 = (2+3*xi-xi.^3)/4;
    N4 = (-1-xi+xi.^2+xi.^3)/4;
    Nu = [N1 le*N2/2 N3 le*N3/2];
    xnume = [xnume;Nx*x_n'];
    unume = [unume;Nu*u_n];
end 

%Exact Solution of curved beam

Ua = -P*R/(2*E*A) + P*R^3/(2*E*I);
Va = -P*pi*R/(4*E*A) - P*pi*R^3/(4*E*I);
Thetaa = -P*R^2/(E*A);

% Plotting
plot(xnume,unume, 'b-', 'LineWidth', 2);
title('Deflection due to point load');
xlabel('Position along the beam (xnume)');
ylabel('Deflection(unume)');
grid on;

disp(size(K))
disp(Freac)
disp(ureduce)
