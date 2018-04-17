clear;clc;clf;

dimension = 3;

data = importdata('data.mat');
%data = importdata('try_data.mat');

X = data(:,1:8);
y = data(:,9);
%X = data(:,1:4);
%y = data(:,5);

[COEFF,SCORE,latent,tsquare, explained, mu] = pca(X);

X = X * COEFF(:,1:dimension);% + repmat(mu(1:dimension), 768, 1);

disp(sum(latent(1:dimension))/sum(latent));

if dimension==2 || dimension==3
    hold on;
    for i=(0:25)
        X0 = X(y==i, :);
        if dimension==2
            plot(X0(:, 1), X0(:, 2), '.', 'color', [rand() rand() rand()]);
        elseif dimension == 3
            plot3(X0(:, 1), X0(:, 2), X0(:, 3), '.', 'color', [rand() rand() rand()]);
        end       
    end
    hold off;
end