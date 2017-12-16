% @author: Kun Zhang, Polytechnique Montreal, 2017-07-04

mFileErrorCode = 100    % Beginning of the m-file

% --- Set general parameters-----------------------------------------------

% Trnsys settings
dt = 0.25*3600; % time step [s]
n = 1/0.25; % 
Nc = 12*n; % control horizon 
Np = 12*n; % prediction horizon
Nsim = 7*24*n+1; % simulation time is 7 days

% Building physical parameters from Trnsys
UA1 = 274.0/3.6; % heat loss coefficient of zone 1 [W/K]
UA2 = 274.7/3.6; 
Cp1 = 13262.057*1000; % building capacitance of zone 1 [J/K]
Cp2 = 4898.682*1000; 
Up = 148.716/3.6; % conductance between 2 zones, [W/K]

% Trnsys inputs, variables of time
msa1 = 0; % ventilation flowrate of zone 1, kg/s
msa2 = 0;
Tsa1 = 20; % ventilation temperature C
Tsa2 = 20;
minf1 = 77.76/3600; % infiltration flowrate of zone 1, kg/s
minf2 = 77.76/3600;
Qig1 = 720/3.6; % internal gains, W
Qig2 = 720/3.6*5; 

% Power price for peak and off-peak periods   
H = cell(Nsim);
for k = 1:n*24                    
    if  ((k >= 5.5*n+1 && k<= 9.5*n+1) || (k >= 16.5*n+1 && k<=20.5*n+1))            
        H{k} = 2;            
    else            
        H{k} = 1;            
    end
end

for k = n*24+1:1:Nsim+Np
    H{k} = H{k-n*24};
end

% Tlb,Tub are the lower and upper band of thermal comfort temperature range
for k = 1:n*24           
    if (k >= 6.5*n+1 && k<=8*n+1 || (k >= 16.5*n+1 && k<=22*n+1))            
        Tlb{k} = 20;
        Tub{k} = 23;            
    else            
        Tlb{k} = 18;
        Tub{k} = 25;
    end        
end
    
for k = n*24+1:1:Nsim+Np
        Tlb{k} = Tlb{k-n*24};
        Tub{k} = Tub{k-n*24};
end

mFileErrorCode = 110    % After setting parameters


%% --- Process Inputs -----------------------------------------------------

Ta = trnInputs(1);
Qsg1 = trnInputs(2);
Qsg2 = trnInputs(3);
Ttrn1 = trnInputs(4);
Ttrn2 = trnInputs(5);

mFileErrorCode = 120    % After processing inputs

% --- First call of the simulation: initial time step (no iterations) -----

if ( (trnInfo(7) == 0) && (trnTime-trnStartTime < 1e-6) )

    % This is the first call (Counter will be incremented later for this 
    % very first call)
    nCall = 0;
    
    % This is the first time step    
    nStep = 1;    
    
    % pre-allocate history of the variables
    nTimeSteps = (trnStopTime-trnStartTime)/trnTimeStep + 1; 
    history.Ta = zeros(nTimeSteps,1); % ambient temperature
    history.Qsg1 = zeros(nTimeSteps,1); % solar gains of zone1
    history.Qsg2 = zeros(nTimeSteps,1); % solar gains of zone2   
    history.Tmat = zeros(nTimeSteps,2); % Matlab instantaneous zone temp
    history.Tmat(1,:) = [20; 20]; % initial value
    history.Tmat_avg = zeros(nTimeSteps,2); % Matlab average zone temp 
    history.Tmat_avg(1,:) = [20; 20]; % initial value
    history.Ttrn = zeros(nTimeSteps,2);
    history.Ttrn(1,:) = [20; 20]; % initial value
    history.x = zeros(nTimeSteps,Nc); % optimal U trajectory
    
    mFileErrorCode = 130    % After initialization
        
end


% Very last call of the simulation (after the user clicks "OK"): Do nothing 

if ( trnInfo(8) == -1 )    
    mFileErrorCode = 1000; 
    
    % draw plot here
end


% if (trnInfo(13) == 1)    
mFileErrorCode = 140;   % Beginning of a post-convergence call
nStep = (trnTime-trnStartTime)/trnTimeStep+1;
    
% return  % Do not update outputs at this call
% end

% --- All iterative calls -------------------------------------------------
Ttrn = history.Ttrn(nStep,:)';
% dry air specific heat in Trnsys for temp lower than 400K, [J/kg.K]
Tak = Ttrn+273.15;
cpa = (0.81764+0.0017855*Tak-0.0000056517*Tak.^2+6.0054*10^-9*Tak.^3)*1000;

% matrices of continuous functions
A = [-(UA1+msa1*cpa(1)+minf1*cpa(1)+Up)/Cp1   Up/Cp1;
         Up/Cp2  -(UA2+msa2*cpa(2)+minf2*cpa(2)+Up)/Cp2];
B = [1/Cp1; 0];
E = [msa1*cpa(1)/Cp1 0 (UA1+minf1*cpa(1))/Cp1 1/Cp1 0 1/Cp1 0; 
     0 msa2*cpa(2)/Cp2 (UA2+minf2*cpa(2))/Cp2  0 1/Cp2 0 1/Cp2];

Ad = expm(A*dt);
Bd = A\(Ad-eye(2))*B;
Ed = A\(Ad-eye(2))*E;
Cd = [1 0];

[Ap,Bp,Ep] = PredictiveMatrix(Ad,Bd,Ed,Cd,Nc,Np);
% perfect forecast of disturbances [Tsa1;Tsa2;Ta;Qsg1;Qsg2;Qig1;Qig2]
load Weather
Tsa1_p = ones(Np,1)*20;
Tsa2_p = ones(Np,1)*20;
Ta_p = Weather(nStep+1:nStep+Np,2);
Qsg1_p = Weather(nStep+1:nStep+Np,3)/3.6;
Qsg2_p = Weather(nStep+1:nStep+Np,3)/3.6;
Qig1_p = Qig1*ones(Np,1);
Qig2_p = Qig2*ones(Np,1);
wp = [Tsa1_p Tsa2_p Ta_p Qsg1_p Qsg2_p Qig1_p Qig2_p]';
Wp = reshape(wp,Np*7,1); % predictive disturbance

Tmat = history.Tmat(nStep,:)';
% Linprog function formulation
Tmin = [Tlb{nStep+1:nStep+Np}]'; % lower bound of output
Tmax = [Tub{nStep+1:nStep+Np}]'; % upper bound of output
f = [H{nStep+1:nStep+Nc}]'; % Power price matrix
Aineq = [-Bp; Bp]; % inequality Ax <= b
bineq = [-Tmin+Ap*Tmat+Ep*Wp; Tmax-Ap*Tmat-Ep*Wp];
Aeq = [];
beq = [];
lb = zeros(Nc,1);
ub = ones(Nc,1)*5160; % maximum power is 5.16kw
options = optimoptions('linprog','Algorithm','interior-point');
[x,fval,~,output,lambda] = linprog(f,Aineq,bineq,Aeq,beq,lb,ub,options);
Qop1 = x(1); % optimal heating rate for zone 1
% calculate average temperature
Ada = 1/dt*inv(A)*(expm(A*dt)-eye(2));
Bda = 1/dt*inv(A)*inv(A)*(expm(A*dt)-eye(2))*B-A\B;
Eda = 1/dt*inv(A)*inv(A)*(expm(A*dt)-eye(2))*E-A\E;
Tmat_avg = Ada*Tmat+Bda*Qop1+Eda*wp(:,1);
% update instantenous temperature
Tmat = Ad*Tmat+Bd*Qop1+Ed*wp(:,1);
PowerPrice = [H{nStep:nStep}]; 
ComfortLb = [Tlb{nStep:nStep}];
ComfortUb = [Tub{nStep:nStep}];

% --- If this is a first call in the time step, increment counter ---

if ( trnInfo(7) == 0 )
    nStep = nStep+1;
end       
    
% --- Get TRNSYS Inputs ---------------------------------------------------

nI = trnInfo(3);   % For bookkeeping
nO = trnInfo(6);   % For bookkeeping

Ta = trnInputs(1);
Qsg1 = trnInputs(2);
Qsg2 = trnInputs(3);
Ttrn1 = trnInputs(4);
Ttrn2 = trnInputs(5);

mFileErrorCode = 150;   % After reading inputs 

history.Ta(nStep) = Ta;
history.Qsg1(nStep) = Qsg1;
history.Qsg2(nStep) = Qsg2;
history.Tmat(nStep,:) = Tmat';
history.Tmat_avg(nStep,:) = Tmat_avg';
history.Ttrn(nStep,:) = [Ttrn1;Ttrn2];
history.x(nStep,:) = x';

% --- Set outputs --------------------------------------------------------

trnOutputs(1) = Tmat_avg(1);
trnOutputs(2) = Tmat_avg(2);
trnOutputs(3) = PowerPrice;
trnOutputs(4) = ComfortLb;
trnOutputs(5) = ComfortUb;
trnOutputs(6) = Qop1*3.6;
trnOutputs(7) = Tmat(1);
trnOutputs(8) = Tmat(2);

% Tell TRNSYS we reached the end of the m-file without errors
mFileErrorCode = 0; 
return  % Do not update outputs at this call
