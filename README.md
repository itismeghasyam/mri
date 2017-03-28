# mri
Compressed sensing in MRI: Study of undersampling patterns, 3D MR image compression using PCA and Optimal Gradient waveforms for a given trajectory

Undersampling_1D.m - Experiment with undersampling pattern constant in 1D (across rows / columns in k-space)

Undersampling_2D.m - Experiment with undersampling pattern constant in 2D (A grid in K-space)

my_cdf.m - Given a vector x with postive elements, spouts out the cdf of the vector by sorting it and normalizing the maximum value to one

ft.m, ift.m - Files provided by Dr.Wong for Fourier transform and Inverse Fourier transform

imshow3Dfull.m, ImageRecon.m - Files provided by Rui, to visualize a 3D MR Image

findnpoints.m - find min number of points in the input vector that contribute to tol% of its energy

image_est.m - Project the input image onto the basis U_basis, at the resolution (number of basis vectors) given by n_basis_vec

pca_compression.m - Code to test the principal component analysis compression technique for a given 3D MR data

