function points_vec = findnpoints(x_vec,tol);

% converting x_vec into a column vector
x_vec = x_vec(:); 

total_energy = x_vec'*x_vec;
x_vec_abs = abs(x_vec);
x_vec_enr = x_vec_abs.^2;

[val pos] = sort(x_vec_enr);

% doing a linear search for min no of points needed to have atleast tol% of
% energy of the signal

for j=1:length(x_vec)
    if sum(val(1:j))> (1-tol)*total_energy
        n_points = j;
        break
    end
end

points_vec = pos(j+1:end);    
    
    
end