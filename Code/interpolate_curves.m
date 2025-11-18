function curve_target = interpolate_curves(aa,a0)

C = cell(20, 1);  
amax = 32;
A = aa';

for k = 1:20
    filename = sprintf('AvsN%d.txt', k);
    C{k} = load(filename);
    C{k}(:,1) = C{k}(:,1) * 4;
end

for k = 1:numel(A)

    indices1 = A(:, 1) <= a0;
    indices2 = A(:, 1) >= a0;
    
    matching_values1= A(indices1, :);
    matching_values2= A(indices2, :);
    
    if a0>=aa(1) && a0<=aa(end)
        index1 = find(A == matching_values1(end));
        index2 = find(A == matching_values2(1));
    elseif a0 < aa(1)
        index1 = 1;
        index2 = 1;
    elseif a0 > aa(end)
        index1 = 20;
        index2 = 20;
    end

end

curve1 = C{index1};
curve2 = C{index2};

n_points = 200; 

a1 = linspace(curve1(1,2),amax,n_points); 
N1 = interp1(curve1(:,2),curve1(:,1),a1,'spline');
a2 = linspace(curve2(1,2),amax,n_points); 
N2 = interp1(curve2(:,2),curve2(:,1),a2,'spline');

a_target = linspace(a0,amax,n_points);
if a0 > aa(end)
    weight = a0/aa(index2);
    N_target = 1/weight *N1;
elseif a0 < aa(1)
    weight = a0/aa(index2);
    N_target = 1/weight *N1;
else
    weight_right = abs((a0 - aa(index1)) / (aa(index2) - aa(index1)));
    weight_left = 1 - weight_right;
    N_target = weight_left *N1 + weight_right*N2;
end

curve_target(:,1) = N_target;
curve_target(:,2) = a_target;

end
