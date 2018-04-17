clear;clf;
n = 3000; %amount of data
m = 2;  %nubmer of attributions
U = 5;
x = fix(n/3);

% begin to generate first cluster
data = zeros(n, m);
for i = 1:m
    data(:, i) =  randn(n, 1)./U;
end 
sum = zeros(n,1);
for i = 1:m
    sum = sum + data(:,i).^2;
end

data = data(find(sum < 1),:);
data1 = data(1:x, : ) * 10;
if m==2,
    hold on;
    plot(data1(:, 1), data1(:, 2), '.');
end

% begin to generate second cluster
data = zeros(n, m);
for i = 1:m
    data(:, i) =  randn(n, 1)./U;
end 
sum = zeros(n,1);
for i = 1:m
    sum = sum + data(:,i).^2;
end

data = data(find(sum < 1),:);
data2 = data(1:x, :) * 10;
for i = 1:fix(x/2),
    p = 0;
    for j = 1:m,
        p = p + data2(i, j).^2;
    end
    p = 10/sqrt(p) + 1;
    for j = 1:m,
        data2(i,j) = data2(i,j) * p;
    end
end

for i = fix(x/2)+1:x,
    p = 0;
    for j = 1:m,
        p = p + data2(i, j).^2;
    end
    p = -10/sqrt(p) - 1;
    for j = 1:m,
        data2(i,j) = data2(i,j) * p;
    end
end

if m==2,
    plot(data2(:, 1), data2(:, 2), '.', 'color', [1 0 0]);
end

% begin to generate third cluster
x1 = n - x*2;
data = zeros(n, m);
for i = 1:m-1
    data(:, i) =  randn(n, 1)./U;
end 
data(:, m) = rand(n, 1)*40;
sum = zeros(n,1);
for i = 1:m-1
    sum = sum + data(:,i).^2;
end

data = data(find(sum < 1),:);
data3 = data(1:x1, :);
p = 25;
data3(:, 1:m-1) = data3(:, 1:m-1) * 5 + p;
data3(:, m) = data3(:, m) - 5;

if m==2,
    plot(data3(:, 1), data3(:, 2), '.', 'color', [0 1 0]);
end

if m==2,
    hold off;
end

data = [data1, zeros(x, 1) ;data2, ones(x, 1) ; data3, ones(x1, 1) * 2];