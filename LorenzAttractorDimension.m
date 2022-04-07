clear all; close all; clc; 
disp('Completion: 0%');

%% Increase rho and record correlation dimension at each step
rhomin = 28; % min, max and number of incriments of rho
rhomax = 50;
numrho = 50; 
rhos = linspace(rhomin,rhomax,numrho);
for r = 1:numrho
    %% Define Symbolic Lorenz Equations
    sigma = 10;
    beta  =8/3;
    
    rho =rhos(r); % Fixed points become repulsors after 24.74

    lorenz = @(t, x)([sigma*(x(2)-x(1)); rho*x(1)-x(2)-x(1)*x(3); x(1)*x(2)-beta*x(3)]);
    
    %% Define Solving conditions
    tspan = [0 100];
    while tspan(2) <= 3
        tspan(2) = input('Input tspan greater than 3: ');
    end
    
    initialConditions = [0,1,0]; % Set starting point
    options = odeset('reltol', 1e-6, 'abstol', 1e-8); % Set error tolerances
    %% Solve System
    [time, soln]=ode45(lorenz, tspan, initialConditions,options ); % solve system of equations
    
    settleTime = 5000;
    points = soln(settleTime:length(soln),:); % Array of solutions after some settling time
    
    %% Increase size of epsilon-ball and record mean & error at each step
    emin = 1; % incriments of size of epsilon-ball in interval 10^emin to 10^emax
    emax = 1.8;
    nume = 1000; 
    epsilons = logspace(emin,emax,nume); % logarithmically spaced incriments in size
    
    C = zeros([1,length(epsilons)]); % pre-load array
    Cerror = zeros([1,length(epsilons)]); % pre-load array
    
    for e=1:length(epsilons)
        epsilon = epsilons(e); % current size of ball
        %% Find number of points within epsilon-ball many times, then average result over multiple points
        sampleNumber = 1000; % number of samples taken for mean
        samples = zeros([1,sampleNumber]); % pre-load array
        for n=1:sampleNumber
            rand = randi(length(points)); % Select one point at random 
            X = points(rand,:);
            
            displacements = zeros([length(points),1]); % pre-load array
            for p=1:length(points)
                % Create array of distances of all points to selected point,should be equal
                % to 0 @ X
                displacements(p) = norm(points(p,:)-X);
            end
            % Find all distances less than or equal to radius 'epsilon'         
            ballElements = find(displacements<=epsilon); %Elements of 'points' that are within ball radius
            samples(n) = length(ballElements);
        end
        
        C(e) = mean(samples); % mean number of elements within epsilon-ball
        Cerror(e) = std(samples)/sqrt(sampleNumber); % standard error in mean
    
        disp(strcat('Completion: ',num2str(100*(r-1)/numrho+(100*e/(numrho*nume))),'%')); % report progress
    end   

    %% Plots and analysis  
    figure 
    errorbar(epsilons,C,Cerror), grid on, hold on;
    set(gca, 'XScale','log', 'YScale','log')
    xlabel('Epsilon [log scale]'), ylabel('C(e) [log scale]')
    title(strcat('Mean point number for increasing size of epsilon-ball, rho=', num2str(rho))); 
    
    % Line of Best fit within scaling region
    lineRange = find(epsilons >= 0.3 & epsilons <= 7.0);
    logx = log(epsilons(lineRange)); % rescale e & C for log
    logy = log(C(lineRange));
    [Const,S] = polyfit(logx, logy, 1); % fitting line
    r2 = 1 - (S.normr/norm(logy - mean(logy)))^2; %r^2 of line of best fit
    
    corDim(r,1) = Const(1); % storing correlation dimension for certain rho
    corDim(r,2) = r2; % with r^2

    plot(epsilons(lineRange), exp(polyval(Const, logx)),'linewidth',2);
    legend('C(e) with standard error in mean',strcat('Least-squares line: m=',num2str(Const(1)),', r^2=',num2str(r2)));
end


figure
plot(rhos,corDim(:,1)), hold on;
xlabel('Rho'), ylabel('C(e)')
title('Change in correlation dimension for increasing rho');

figure 
figure(1)
comet3(soln(:,1), soln(:,2), soln(:,3))
xlabel('x');ylabel('y');zlabel('z')
title('Time Evolution of solutions to Lorenz Equations');