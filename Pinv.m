% The Moore-Penrose pseudo inverse for the matrix (M - lambda * I)
%
% Authors: Jin Wu, Yuxiang Sun, Miaomiao Wang, Ming Liu
% e-mail: jin_wu_uestc@hotmail.com
% website: www.jinwu.science
%          www.ram-lab.com
%
% Reference: Wu, J., Sun, Y., Wang, M., Liu, M. (2019) 
%               Hand-eye Calibration: 4D Procrustes Analysis Approach,
%               Submitted to IEEE Trans. Instrum. Meas.


function B = Pinv(A)
[u, s, v] = svd(A);
B = v * diag([1 / s(1, 1), 1 / s(2, 2), 1 / s(3, 3), 0]) * u';
end