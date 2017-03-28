function I = image_est(img_in, U_basis,n_basis_vec)

img_vec = img_in(:);
ci_vec = zeros(1, n_basis_vec);


for j = 1:n_basis_vec
    ci_vec(j) = img_vec'*U_basis(:,j)/(norm(U_basis(:,j)))^2;
end

img_est = zeros(size(img_vec));

for j = 1:n_basis_vec
    img_est = img_est + ci_vec(j)*U_basis(:,j);
end

I = reshape(img_est, [size(img_in,1) size(img_in,2)]);




end