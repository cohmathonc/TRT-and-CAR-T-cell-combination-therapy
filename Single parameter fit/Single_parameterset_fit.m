close all
clear all
tic
colormap hsv
global rho nParams lambda_p gamma TV0 time n_cohort cohorts p_current_set inj_act tCAR N_C0 tTRT cohort_count N_C_transf

cur_dir = pwd;
cd ..
data_dir = pwd;
cohort_labels = {'TRT D7','TRT D7 + CAR-T D18','TRT D7 + CAR-T D25', 'TRT D7 + CAR-T D32','CAR-T D7'}
colormap hsv
gamma = 0.693/(60/(24*60));                  % Lee-Catcheside decay factor assuming an interaction window of 16 mins for double strand break. Mazeron et al book chapter and dale book chapter ~ 1 hr
N_C0 = [0,1e6,1e6,1e6,1e6];     %initial number of CAR T cells for each cohort
tCAR = [7,18,25,32,7];        % CAR-T Rx dates
tTRT = [7,7,7,7,0];
N_C_transf = zeros(size(N_C0));

cohorts = {'Ac_200nCi','Ac_200nCi+CAR-T D18','Ac_200nCi+CAR-T D25','Ac_200nCi+CAR-T D32','CAR-T only'};
inj_act = [0.2,0.2,0.2,0.2,0]; % uCi
rho_0 = 0.27;
eta_0 = 4;               % Fixed from first cohort fit
k2_0 = 1e-12;             % /day/cell     % in base case - 1e-11
alpha_T_0 = 1.5;            % /Gy
alpha_C_0 = 1.2*alpha_T_0;      % /Gy
k1_0 = 1e-6;              % /day/cell
theta_0 = 0.02;             % /day
kcl_0 = .6;                % /day
rho = rho_0;
DF_d7 = exp(rho_0*7);       % decay factor to day 7


n_cohort = 5;
lambda_p = 0.693/9.92;                    % decay constant in /sec assuming t1/2 = 9.92d

matrix_count = 0;
mouse_count = 0;

%% READ IN THE TUMOR VOLUME AND TIME DATA INTO ARRAYS TV AND TIME_MATRIX AND INITIALIZE TUMOR VOLUMES FOR EACH MOUSE
cd(data_dir);
vol_data = Load_data('all_cohorts');
cd(cur_dir);

for cohort_count = 1:n_cohort
    time = vol_data{cohort_count}(:,1)';   TV = vol_data{cohort_count}(:,2:end)';
    n_mouse = size(TV,1);
    for mouseNum = 1:n_mouse
        valid_ind = min(find(isnan(TV(mouseNum,:))))-1;         % If a mouse dies data further from that time is NaN - ignore the NaNs (valid_ind+1:end)
        if isempty(valid_ind)
            valid_ind = length(TV(mouseNum,:));
        end
       
    % time_matrix will store time stamps of BLI data for each mouse
    % formatted in a columnwise format. - Time - Mouse # - Cohort #
    TV_matrix(matrix_count+1:matrix_count+valid_ind) = TV(mouseNum,1:valid_ind);
    time_matrix(matrix_count+1:matrix_count+valid_ind,1) = time(1:valid_ind);           % time point of data post tumor innoculation in hours
    time_matrix(matrix_count+1:matrix_count+valid_ind,2) = mouse_count+mouseNum;        % label the mouse number - used for parameter optimization
    time_matrix(matrix_count+1:matrix_count+valid_ind,3) = cohort_count;                % cohort number - used for parameter optimization
    matrix_count = matrix_count + valid_ind;
    TV0(mouse_count+mouseNum) = TV(mouseNum,1);
    end
    mouse_count = mouse_count +mouseNum;
end
NMice = mouse_count;    clear mouse_count;
BLI2Nt = 3.72;          % from control cohort
TV_matrix = BLI2Nt*TV_matrix;
TV0 = BLI2Nt*TV0;


%% Initialize the parameters
% p0 is the initial parameter set that is input to the optimization loop

nGlobalParams = 8;      % eta0 and k2 and alpha_C;
nParams = 0;            % (# of params per mouse) - alpha_T, k1, theta, k_cl

p0(1) = eta_0;  lb(1) = 0.45;  ub(1) = 10;                  % eta - upper limit depends on the biological clearance. Lower limit = lambda_physical
p0(2) = k2_0;    lb(2) = k2_0/1000;  ub(2) = k2_0*1000;      %k2 for CAR-T cells (prolif/exhaustion rate)
p0(3) = alpha_C_0;  lb(3) = .05;  ub(3) = 5;                %alpha_C for all mice is same - same CAR-T cells
    
p0(4) = alpha_T_0;  lb(4) = 1;  ub(4) = 2;            %     alpha_T
p0(5) = k1_0;  lb(5) = 1e-8;  ub(5) = k1_0*1000;        % k1
p0(6) = theta_0;  lb(6) = 0.0;  ub(6) = .1;              % theta
p0(7) = kcl_0;  lb(7) = 0.05;  ub(7) = 1;                % clearance constant              % theta
p0(8) = rho_0;  lb(8) = 0.1;    ub(8) = 0.35;       %rho

%%  RUN OPTIMIZATION FOR THE PARAMETER SET
options = optimoptions('lsqcurvefit','MaxFunctionEvaluations',10000,'MaxIterations',30000);%,'FunctionTolerance',1e-15,'StepTolerance',1e-10);
TV_matrix = log(TV_matrix);
[p,resnorm] = lsqcurvefit(@RIT,p0,time_matrix,TV_matrix,lb,ub,options)
TV_matrix = exp(TV_matrix);
colr = lines(max(time_matrix(:,2)));
k_cohort = NaN(10,n_cohort);    alpha_cohort = k_cohort;   k_cl_cohort = k_cohort;

hold on
    figure('Color','w','Position',[0,0,1500,850]);
%% PLOT THE RESULTS
for cohort_count = 1:n_cohort
    mouseNum_cohort = 1;
    subplot(2,3,cohort_count);
    ind_first = find(time_matrix(:,3) == cohort_count,1,'first');
    ind_last = find(time_matrix(:,3) == cohort_count,1,'last');

    for mouseNum = time_matrix(ind_first,2):time_matrix(ind_last,2)
        ind = find(time_matrix(:,2) == mouseNum);
        time = time_matrix(ind,1);
        
        p_current_set = p;%[p(1:3),p((mouseNum-1)*nParams+nGlobalParams+1:(mouseNum-1)*nParams+nGlobalParams+4)]
        
         ind_CART = find(time>tCAR(cohort_count),1);
         if ~N_C0(cohort_count)
             soln = ode45(@DifEq, time, [TV0(mouseNum);0;0],odeset('MaxStep',0.1));                                              % calculate the solution                                   
             t = time(1):time(end);
             y = deval(soln,t);
             y = sum(y(1:2,:));
         else
                soln_TRT = ode23(@DifEq, time(1:ind_CART), [TV0(mouseNum);0;0]);
                vol_TRT = deval(soln_TRT,time(1):tCAR(cohort_count));
                y0_CART = vol_TRT(:,end)+[0;0;N_C0(cohort_count)];
                %add the CART at this point
                soln_CART_TRT = ode23(@DifEq, [tCAR(cohort_count);time(ind_CART:end)], y0_CART);
                vol_TRT_CART= deval(soln_CART_TRT,tCAR(cohort_count):time(end));
                y = [sum(vol_TRT(1:2,1:end-1)),sum(vol_TRT_CART(1:2,:))];
                t = time(1):time(end);
         end

        semilogy(time,TV_matrix(ind),'s','color',colr(mouseNum,:),'MarkerSize',3,'MarkerFaceColor',colr(mouseNum,:));  % plot on log axis
        hold on
        semilogy(t,y,'color',colr(mouseNum,:),'LineWidth',4);                % plot on log axis
        ind_pfs = max(find(~(y>TV0(mouseNum))))+1;
        pfs(cohort_count,mouseNum_cohort) = t(ind_pfs) - min(tTRT(cohort_count),tCAR(cohort_count));
     
        [~,ind_min] = min(y);
        tmin(cohort_count,mouseNum_cohort) = t(ind_min)- min(tTRT(cohort_count),tCAR(cohort_count));
%        PFS(mouseNum) = t(find((y>TV0(mouseNum)).*(t>tmin(mouseNum)).*(t>tCAR(cohort_count)),1));
%        OS(mouseNum) = t(find(y>1e11,1));
%          
%        semilogy(time(1):time(end),deval(soln,time(1):time(end)),'--','color',colr(mouseNum,:),'LineWidth',2);                % plot on log axis
%        min_tv(mouseNum) = min(sum(soln.y));                                         % minimum tumor burden
%        min_tv_t(mouseNum) = soln.x(find(sum(soln.y)==min_tv(mouseNum)))-time(1)                      % time @ minimum tumor burden
        ind_min_tv = max(find(diff(sum(soln.y))<0));
        if ~isempty(ind_min_tv)
            inflection_t(cohort_count,mouseNum_cohort) = soln.x(ind_min_tv)-time(1);
        end
        mouseNum_cohort = mouseNum_cohort+1;
    end
   
    box off
    title(cohort_labels(cohort_count))
    xlabel('Time post innoculation (days)');        ylabel('Tumor burden (cells)');       set(gca,'FontSize',18,'LineWidth',1);   axis([0,100,1e5,3e11]);    
    
   
end

pfs(5,end) = NaN;       tmin(5,end) = NaN;
subplot(2,3,6)
h = boxplot(pfs','labels',cohort_labels);   hold on
h1 = boxplot(tmin','labels',cohort_labels,'Color','b');
set(h,{'linew'},{2});   box off
set(h1,{'linew'},{2});   box off
ylim([0,60]);
ylabel('Survival parameter (days)')
title('Survival');
set(gca,'FontSize',18,'LineWidth',2);


%% FUNCTION TO BE OPTIMIZED - CONTAINS TWO LOOPS RUNNING OVER MICE AND COHORTS TO CALCULATE THE SOLUTION FOR EACH ITERATION OF PARAMETER SET
function S = RIT(p, time_matrix, ~)
    global lambda_p TV0 n_cohort nParams inj_act rho tCAR N_C0 tTRT 
    
    y0 = TV0;

    for cohort_count = 1:n_cohort
      
        %tdata is initialized every time a new cohort is run
        tdata = time_matrix(find(time_matrix(:,3) == cohort_count),:);                  % use time data for only the mice in the current cohort
        firstMouseNum = tdata(1,2,1);   lastMouseNum = tdata(end,2,1);                  % mice in the current cohort
        
      % For each mouse number, calculate the solution of tumor volume curve given
      % the current parameters for that mouse. Solution goes into FitData.
      for datacount = firstMouseNum:lastMouseNum
            ind = find((tdata(:,2) == datacount)&(tdata(:,3) == cohort_count));         % indices of time for current mouse
            ind_CART = find(tdata(ind,1)>tCAR(cohort_count),1);                         %first index after CAR-T treatment start
            
            if isempty(ind_CART)
                FitData{datacount} = ode23(@DifEq, tdata(ind,1), [y0(datacount);0;0]);
                S{datacount} = deval(FitData{datacount},tdata(ind,1));                      % Calculate the tumor volume at the defined time points of input data using deval function and the generated fit
            else
                %use data up to ind_CART data point for the TRT phase since
                %the data in between might be left without fitting and the
                %initial condition for tCART phase will not be there
                FitData_0{datacount} = ode23(@DifEq, tdata(ind(1:ind_CART),1), [y0(datacount);0;0]);
                S_0{datacount} = deval(FitData_0{datacount},[tdata(ind(1:ind_CART-1),1);tCAR(cohort_count)]);
                y0_CART = S_0{datacount}(:,end)+[0;0;N_C0(cohort_count)];

                %add the CART at this point
                FitData_1{datacount} = ode23(@DifEq, [tCAR(cohort_count);tdata(ind(ind_CART:end),1)], y0_CART);
                S_1{datacount} = deval(FitData_1{datacount},[tCAR(cohort_count);tdata(ind(ind_CART:end),1)]);
                S{datacount} = [S_0{datacount}(:,1:end-1),S_1{datacount}(:,2:end)];
                clear S_0 S_1
            end

       end
        
    end
    % S contains the evaluated tumor volumes for each mouse at the given input time points
    S = cell2mat(S);    %max(S(3,:))
    S  = sum(S(1:2,:));
    S = log(S);
    %% FUNCTION FOR CALCULATING THE DIFFERENTIAL CHANGES IN EACH LOOP FOR OPTIMIZATION
% simulate till the time of injection with TRT
% and then simulate with TRT + CART after reinitializing y
    function dydt = DifEq(t, y)   

        N_T = y(1);
        N_R = y(2);
        N_C= y(3);

        R0 = p(1)*inj_act(cohort_count);
        k2 = p(2);
        alpha_C = p(3);
        alpha_T = p((datacount-1)*nParams+4);
        k1 = p((datacount-1)*nParams+5);
        theta = p((datacount-1)*nParams+6);
        k_cl = p((datacount-1)*nParams+7);
        rho = p((datacount-1)*nParams+8);
        
        
        kRx_T = alpha_T*R0*exp(-lambda_p*(t-tTRT(cohort_count)))*(t>tTRT(cohort_count));
        kRx_C = alpha_C*R0*exp(-lambda_p*(t-tTRT(cohort_count)))*(t>tTRT(cohort_count));
        
        dydt = [rho*N_T - kRx_T*N_T - k1*N_T*N_C;...        % N_T
                kRx_T*N_T - k1*N_R*N_C - k_cl*N_R;...            % N_R
                k2*(N_T + N_R)*N_C - kRx_C*N_C - theta*N_C];         %CAR-T cells @ tumor site
       
    end
end

%% FUNCTION FOR CALCULATING THE DIFFERENTIAL CHANGES IN EACH LOOP FOR OPTIMIZATION
% This function is used for input to deval.

function dydt = DifEq(t, y)   
    global lambda_p time p_current_set tTRT cohort_count inj_act rho
    
        N_T = y(1);
        N_R = y(2);
        N_C= y(3);

        R0 = p_current_set(1)*inj_act(cohort_count);
        k2 = p_current_set(2);
        alpha_C = p_current_set(3);
        alpha_T = p_current_set(4);
        k1 = p_current_set(5);
        theta = p_current_set(6);
        k_cl = p_current_set(7);
        rho = p_current_set(8);

        kRx_T = alpha_T*R0*exp(-lambda_p*(t-tTRT(cohort_count)))*(t>tTRT(cohort_count));
        kRx_C = alpha_C*R0*exp(-lambda_p*(t-tTRT(cohort_count)))*(t>tTRT(cohort_count));
        
        dydt = [rho*N_T - kRx_T*N_T - k1*N_T*N_C;...        % N_T
            kRx_T*N_T - k1*N_R*N_C - k_cl*N_R;...            % N_R
            k2*(N_T + N_R)*N_C - kRx_C*N_C - theta*N_C];        %CAR-T cells @ tumor site
      
    T = t-time(1);
    
 end


