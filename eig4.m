% The generalized analytical eigen-decomposition method for symmetric 4x4 matrix
%
% Author: Jin Wu 
% e-mail: jin_wu_uestc@hotmail.com
% website: www.jinwu.science
%
% Reference: Wu, J. et al. (2019) Fast Symbolic 3D Registration Solution.
%                 IEEE Trans. Auto. Sci. Eng. 

function [v, lambda] = eig4(A)
a11 = A(1, 1); a12 = A(1, 2); a13 = A(1, 3); a14 = A(1, 4);
a22 = A(2, 2); a23 = A(2, 3); a24 = A(2, 4);
a33 = A(3, 3); a34 = A(3, 4);
a44 = A(4, 4);

a12_2 = a12 * a12; a13_2 = a13 * a13; a14_2 = a14 * a14;
a23_2 = a23 * a23; a24_2 = a24 * a24;
a34_2 = a34 * a34;

tau1 = -a11 - a22 - a33 - a44;
tau2 = -a12_2 - a13_2 - a14_2 + a11*a22 - a23_2 - a24_2 + a11*a33 + a22*a33 - a34_2 + a11*a44 + a22*a44 + a33*a44;
tau3 = a13_2*a22 + a14_2*a22 - 2*a12*a13*a23 + a11*a23_2 - ...
      2*a12*a14*a24 + a11*a24_2 + a12_2*a33 + a14_2*a33 - ...
      a11*a22*a33 + a24_2*a33 - 2*a13*a14*a34 - 2*a23*a24*a34 + ...
      a11*a34_2 + a22*a34_2 + a12_2*a44 + a13_2*a44 - ...
      a11*a22*a44 + a23_2*a44 - a11*a33*a44 - a22*a33*a44;
tau4 = a14_2*a23_2 - 2*a13*a14*a23*a24 + a13_2*a24_2 - ...
   a14_2*a22*a33 + 2*a12*a14*a24*a33 - a11*a24_2*a33 + 2*a13*a14*a22*a34 - ...
   2*a12*a14*a23*a34 - 2*a12*a13*a24*a34 + 2*a11*a23*a24*a34 + ...
   a12_2*a34_2 - a11*a22*a34_2 - a13_2*a22*a44 + ...
   2*a12*a13*a23*a44 - a11*a23_2*a44 - a12_2*a33*a44 + a11*a22*a33*a44;

tau1_2 = tau1 * tau1;
a = tau2 - 3 / 8 * tau1_2;
b = tau3 - tau1 * tau2 / 2 + tau1_2 * tau1 / 8;
c = tau4 - tau1 * tau3 / 4 + tau1_2 * tau2 / 16 - 3 * tau1_2 * tau1_2 / 256;

T0 = 2 * a^3 + 27 * b^2 - 72 * a * c;
theta = atan2(sqrt(4 * (a * a + 12 * c)^3 - T0 * T0), T0);
aT1 = 1.259921049894873 * sqrt(a * a + 12 * c) * cos(theta / 3);
T2 = sqrt( - 4 * a + 3.174802103936399 * aT1);
    
lambda = - 0.204124145231932 * (T2 + sqrt( - T2 * T2 - 12 * a + 29.393876913398135 * b / T2)) - tau1 / 4;

G11 = A(1, 1) - lambda; G12 = A(1, 2); G13 = A(1, 3); G14 = A(1, 4);
G22 = A(2, 2) - lambda; G23 = A(2, 3); G24 = A(2, 4);
G33 = A(3, 3) - lambda; G34 = A(3, 4);

q = [
		       G14 * G23 * G23 - G13 * G23 * G24 - G14 * G22 * G33 + G12 * G24 * G33 + G13 * G22 * G34 - G12 * G23 * G34;
                       G13 * G13 * G24 + G12 * G14 * G33 - G11 * G24 * G33 + G11 * G23 * G34 - G13 * G14 * G23 - G13 * G12 * G34;
	               G13 * G14 * G22 - G12 * G14 * G23 - G12 * G13 * G24 + G11 * G23 * G24 + G12 * G12 * G34 - G11 * G22 * G34;
	               - (G13 * G13 * G22 - 2 * G12 * G13 * G23 + G11 * G23 * G23 + G12 * G12 * G33 - G11 * G22 * G33)];
v = q ./ norm(q);
end