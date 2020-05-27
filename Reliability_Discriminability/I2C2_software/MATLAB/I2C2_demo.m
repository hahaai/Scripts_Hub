%Example usage of I2C2 functions

%import the data into MATLAB
y = importdata('vn10_t5.txt');

%compute the I2C2 coefficient
[vn_lambda,~,~,vn_demean_y] = I2C2_compute(y, 'I',21, 'J',2, 'demean', true);

%compute the confidence intervals
[vn_ci] = I2C2_mcCI(y,'i',21,'j',2)

%compute the null distribution
vn_NullDist=  I2C2_mcNulldist( vn_demean_y, 'I',21, 'J', 2, 'R', 100, 'rseed', 1, 'ncores', 1, 'demean',false )