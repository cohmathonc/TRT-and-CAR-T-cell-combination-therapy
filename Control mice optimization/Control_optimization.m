close all
clear all

global TV0 time n_cohort n_control p_current_set

colormap hsv

% BLI data for control mice tumor growth

input_data = [7	8.99e+006	1.28e+007	9.74e+006	1.08e+007	7.15e+006	7.67e+006	7.57e+006	3.34e+006; ...
14	9.12e+007	1e+008	7.29e+007	5.37e+007	7.52e+007	5.63e+007	8.24e+007	6.67e+006; ...
21	4.76e+007	1.07e+009	7.5e+008	5.09e+008	3.82e+008	5.71e+008	6.87e+008	9.08e+008; ...
35	4.48e+010	4.61e+010	3.48e+010	2.31e+010	2.44e+010	2.79e+010	2.77e+010	3.07e+010; ...
42	5.62e+009	1.91e+011	1.86e+011	1.93e+011	1.21e+011	9.79e+010	1.52e+011	1.56e+011];	

t = input_data(:,1);
BLI = input_data(:,2:end);  % ignoring mouse 1 as the last data point seems to have lower BLI
BLI(1,end) = NaN;


Control_vol = BLI';
Control_time = t;


n_cohort = 1;
matrix_count = 0;               %Keeps track of the number of valid data points that have time and mouse volume
mouse_count = 0;                %Keeps track of the mouse number

% For each cohort, read the tumor volume and time point data in the TV and
% time matrix. Find the number of valid data points (goes in matrix_count)

%% READ IN THE TUMOR VOLUME AND TIME DATA INTO ARRAYS TV AND TIME_MATRIX AND INITIALIZE TUMOR VOLUMES FOR EACH MOUSE

for cohort_count = 1:n_cohort
    n_mouse = size(Control_vol,1);
    for mouseNum = 1:n_mouse
        valid_ind = min(find(isnan(Control_vol(mouseNum,:))))-1;         % If a mouse dies data further from that time is NaN - ignore the NaNs (valid_ind+1:end)
        if isempty(valid_ind)
            valid_ind = length(Control_vol(mouseNum,:));
        end
        
    TV_matrix(matrix_count+1:matrix_count+valid_ind) = Control_vol(mouseNum,1:valid_ind);
    time_matrix(matrix_count+1:matrix_count+valid_ind,1) = Control_time(1:valid_ind);           % time point of data post tumor innoculation in days
    time_matrix(matrix_count+1:matrix_count+valid_ind,2) = mouse_count+mouseNum;        % label the mouse number - used for parameter optimization
    time_matrix(matrix_count+1:matrix_count+valid_ind,3) = cohort_count;                % cohort number - used for parameter optimization
    matrix_count = matrix_count + valid_ind;
    TV0(mouse_count+mouseNum) = Control_vol(mouseNum,1);
    
    end
    mouse_count = mouse_count +mouseNum;
end
n_control = max(time_matrix(find(time_matrix(:,3) == 1),2));        %total number of mice

%% INITIALIZE PARAMETER SET AND THE BOUNDS FOR THE PARAMETERS
%p0 is the initial parameter set for proliferation rate that goes into the lsq fitting algorithm. 

p0(1:n_control) = 0.35;   lb(1:n_control) = 0;  ub(1:n_control) = 1;                    % Proliferation rate rho for control

%%  RUN OPTIMIZATION FOR THE PARAMETER SET
options = optimoptions('lsqcurvefit','MaxFunctionEvaluations',10000,'MaxIterations',30000);
TV_matrix = log(TV_matrix);
[p,resnorm] = lsqcurvefit(@RIT,p0,time_matrix,TV_matrix',lb',ub',options)
r2 = 1-resnorm/sum((TV_matrix'-mean(TV_matrix)).^2);
TV_matrix = exp(TV_matrix);

colr = lines(max(time_matrix(:,2)));            %line colors for plotting
k_cohort = NaN(10,n_cohort);                    % Save the profliferation rate in a matrix
    
%% PLOT THE DATA AND CALCULATED SOLUTION
for cohort_count = 1:n_cohort
    figure('Color','w');
    ind_first = find(time_matrix(:,3) == cohort_count,1,'first');
    ind_last = find(time_matrix(:,3) == cohort_count,1,'last');

    for mouseNum = time_matrix(ind_first,2):time_matrix(ind_last,2)
        ind = find(time_matrix(:,2) == mouseNum);
        time = time_matrix(ind,1);
        p_current_set = p(mouseNum);                                  %for control
        soln = ode45(@DifEq, time, TV0(mouseNum));                                              % calculate the solution                                   

        semilogy(time,TV_matrix(ind),'o','color',colr(mouseNum,:),'MarkerSize',5,'MarkerFaceColor',colr(mouseNum,:));  % plot on log axis
        hold on
        semilogy(time(1):time(end),deval(soln,time(1):time(end)),'-','color',colr(mouseNum,:),'LineWidth',2);                % plot on log axis
    end
    box off
    xlabel('Time post inoculation (days)');        ylabel('Tumor burden (rad)');       set(gca,'FontSize',15,'LineWidth',1);   axis([0,80,1e5,3e11]);    
    
end
dec_fac = exp(mean(p)*7);       %factor accounting for 7 days of exponential growth
BLId0 = nanmean(TV0)/dec_fac;
N0 = 5e6;
BLI2Nt = N0/BLId0;

%% FUNCTION TO BE OPTIMIZED - CONTAINS TWO LOOPS RUNNING OVER MICE AND COHORTS TO CALCULATE THE SOLUTION FOR EACH ITERATION OF PARAMETER SET
function S = RIT(p, time_matrix, ~)
    global  TV0 n_cohort datacount

    y0 = TV0;

    % For each mouse number, calculate the solution of tumor volume curve given
    % the current parameters for that mouse. Solution goes into FitData.
    for cohort_count = 1:n_cohort
    
        tdata = time_matrix(find(time_matrix(:,3) == cohort_count),:);                  % use time data for only the mice in the current cohort
        firstMouseNum = tdata(1,2,1);   lastMouseNum = tdata(end,2,1);

        for datacount = firstMouseNum:lastMouseNum
            ind = find((tdata(:,2) == datacount)&(tdata(:,3) == cohort_count));
            FitData{datacount} = ode45(@DifEq, tdata(ind,1), y0(datacount));
            S{datacount} = deval(FitData{datacount},tdata(ind,1));                      % Calculate the tumor volume at the defined time points of input data using deval function and the generated fit
        end

    end
    % S contains the evaluated tumor volumes for each mouse at the given input time points
    S = cell2mat(S);
    S = reshape(S',size(time_matrix,1),1);
    S = log(S);
    %% FUNCTION FOR CALCULATING THE DIFFERENTIAL CHANGES IN EACH LOOP FOR OPTIMIZATION

    function dydt = DifEq(t, y)   
        rho = p(datacount);
        dydt = rho*y;                 % exponential function
    end
end

%% FUNCTION FOR CALCULATING THE DIFFERENTIAL CHANGES
% This function is used for input to deval. Check if earlier loop can be
% used.

function dydt = DifEq(t, y)   
global p_current_set
    rho = p_current_set;
    dydt = rho*y;
end
