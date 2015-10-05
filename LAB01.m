%load data
y = load('ml1y.dat');
x = load('ml1x.dat');

%Plot
figure, plot(x,y,'rh', 'markerfacecolor', 'auto');
xlabel('Age (yrs.)');
ylabel('Height (m)');

m = length(x);
x = [ones(m, 1), x];

%Procedure 1
    alpha = 0.07;
    theta = zeros(2,1); 

    %1st iteration
    theta = theta - alpha .* ((1/m).* x' * ((x * theta) - y))
        %h = (x * theta)
        %gradient = ((1/m).* x' * (h - y))
        %theta = theta - alpha .* gradient

%Procedure 2
    %Continue running remaining iteration
    for i = 1:1499
               % Dimensions:      2x50 (50x2 * 2x1  - 50x1) = 2x1 
               theta = theta - alpha .* ((1/m).* x' * ((x * theta) - y));
    end
    %Final Theta value
    theta
    
     hold on; %Plot new data without clearing old plot
     plot(x(:,2), x*theta, 'g-'); %x is a matrix with 2 columns, 2nd column containing the time info
     legend('Training Data','Linear Regression')

%Procedure 3
    %Prediction on age = 3.5
     %printf('Height of 3.5 year old boy:')
     height1 = [1, 3.5]*theta
     plot(x, height1*ones(size(y)),'b:');
     plot(3.5*ones(size(y)), y,'b:')

     %Prediction on age = 7
     %printf('Height of 7 year old boy:')
     height2 = [1, 7]*theta
     plot(x, height2*ones(size(y)),'b:');
     plot(7*ones(size(y)), y,'b:')
%Contour Plot
    % Calculate J matrix
    
    % Grid over which we will calculate J
    theta0_vals = linspace(-3, 3, 50);
    theta1_vals = linspace(-1, 1, 50);
    
    % initialize J_vals to a matrix of 0's
    J_vals = zeros(length(theta0_vals), length(theta1_vals));

    for i = 1:length(theta0_vals)
          for j = 1:length(theta1_vals)
          t = [theta0_vals(i); theta1_vals(j)];
          %  50x2 * 2x50 -50x1  50x2 * 2x50 - 50x1
          J_vals(i,j) = (1/(2*m)) .* (x * t - y)' * (x * t - y); %Code for Jtheta formula
          end
    end

    %Plot the surface plot
    J_vals = J_vals'; 
    
    figure, contour(theta0_vals, theta1_vals, J_vals, logspace(-1.5,1.5,12))
    xlabel('\theta_0');
    ylabel('\theta_1');
    
    
    figure;
    surf(theta0_vals,theta1_vals,J_vals)
    xlabel('\theta_0');
    ylabel('\theta_1');
    

