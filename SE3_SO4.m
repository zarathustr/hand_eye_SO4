% The mapping from SE(3) to SO(4)
% 
% Authors: Jin Wu, Yuxiang Sun, Miaomiao Wang, Ming Liu
% e-mail: jin_wu_uestc@hotmail.com
% website: www.jinwu.science
%          www.ram-lab.com
%
% Reference: Wu, J., Sun, Y., Wang, M., Liu, M. (2019) 
%               Hand-eye Calibration: 4D Procrustes Analysis Approach,
%               Submitted to IEEE Trans. Instrum. Meas.


function SA = SE3_SO4(A, dd)
RA = A(1 : 3, 1 : 3);
tA = A(1 : 3, 4);
SA = [RA, tA / dd;
      - tA' * RA / dd, 1];
end