<<<<<<< HEAD
function H = HessFcn(x,lambda,ldata)
global Jac
meas_mids_std = [];
Sigma = [];
for i = 1:length(ldata.emus)
    meas_mid_std = ldata.mids_std{i};
    meas_mid_std = meas_mid_std(~isnan(meas_mid_std));
    Sigma = [Sigma;meas_mid_std];
end
H = Jac'*diag(Sigma)*Jac;
=======
function H = HessFcn(x,lambda,ldata)
global Jac
meas_mids_std = [];
Sigma = [];
for i = 1:length(ldata.emus)
    meas_mid_std = ldata.mids_std{i};
    meas_mid_std = meas_mid_std(~isnan(meas_mid_std));
    Sigma = [Sigma;meas_mid_std];
end
H = Jac'*diag(Sigma)*Jac;
>>>>>>> 2024247f7ad30fed5e7c393c7eddc78b1a8690d7
end