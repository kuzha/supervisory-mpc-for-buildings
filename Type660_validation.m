% This model is to simulate a 2-zone building with. The Trnsys model uses 
% 2 types of Type 660.
% @author: Kun Zhang, Polytechnique Montreal, 2017-07-31

mFileErrorCode = 100    % Beginning of the m-file

% --- Set general parameters-----------------------------------------------
% Trnsys settings
dt = 0.25*3600; % time step [s]

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
Qig2 = 720/3.6; 

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
    history.Tmat(1:2,:) = [20 20; 20 20]; % initial value
    history.Tmat_avg = zeros(nTimeSteps,2); % Matlab average zone temp 
    history.Tmat_avg(1,:) = [20; 20]; % initial value
    history.Ttrn = zeros(nTimeSteps,2);
    history.Ttrn(1,:) = [20; 20]; % initial value
   
    mFileErrorCode = 130    % After initialization
        
end


% Very last call of the simulation (after the user clicks "OK"): Do nothing 

if ( trnInfo(8) == -1 )    
    mFileErrorCode = 1000; 
    
    % draw plot here
end


% --- Post convergence calls: store values --------------------------------

if (trnInfo(13) == 1)
    
    mFileErrorCode = 140;   % Beginning of a post-convergence call 
    nStep = (trnTime-trnStartTime)/trnTimeStep+1;    
    
    history.Ta(nStep)   = Ta;  
    history.Qsg1(nStep)   = Qsg1;
    history.Qsg2(nStep)   = Qsg2;    
    history.Tmat(nStep,:)   = Tmat';
    history.Tmat_avg(nStep,:)   = Tmat_avg';
    history.Ttrn(nStep,:) = [Ttrn1;Ttrn2];
    
    % Tell TRNSYS that we reached the end of the m-file without errors
    mFileErrorCode = 0; 
    return  % Do not update outputs at this call

end

% --- All iterative calls -------------------------------------------------

% --- If this is a first call in the time step, increment counter ---

if ( trnInfo(7) == 0 )
    nStep = nStep+1;
end       
    
% --- Get TRNSYS Inputs ---------------------------------------------------

nI = trnInfo(3);     % For bookkeeping
nO = trnInfo(6);   % For bookkeeping

Ta = trnInputs(1);
Qsg1 = trnInputs(2);
Qsg2 = trnInputs(3);
Ttrn1 = trnInputs(4);
Ttrn2 = trnInputs(5);

mFileErrorCode = 150;   % After reading inputs 

% Main function

% dry air specific heat in Trnsys for temp lower than 400K, [J/kg.K]
Tmat = history.Tmat(nStep-1,:)';
Tak = Tmat+273.15;
cpa = (0.81764+0.0017855*Tak-0.0000056517*Tak.^2+6.0054*10^-9*Tak.^3)*1000;

% matrices of continuous functions
A = [-(UA1+msa1*cpa(1)+minf1*cpa(1)+Up)/Cp1   Up/Cp1;
         Up/Cp2  -(UA2+msa2*cpa(2)+minf2*cpa(2)+Up)/Cp2];
B = [1/Cp1 0; 0 1/Cp2];
E = [msa1*cpa(1)/Cp1 0 (UA1+minf1*cpa(1))/Cp1 1/Cp1 0 1/Cp1 0; 
     0 msa2*cpa(2)/Cp2 (UA2+minf2*cpa(2))/Cp2  0 1/Cp2 0 1/Cp2];
W = [Tsa1;Tsa2;Ta;Qsg1/3.6;Qsg2/3.6;Qig1;Qig2]; % uncontrolled disturbance

% calculate average output value
Ada = 1/dt*inv(A)*(expm(A*dt)-eye(2));
Eda = 1/dt*inv(A)*inv(A)*(expm(A*dt)-eye(2))*E-A\E;
Tmat_avg = Ada*Tmat+Eda*W;

% calculate instantaneous output value
Ad = expm(A*dt);
Ed = A\(Ad-eye(2))*E;
Tmat = Ad*Tmat+Ed*W;

T1_Matlab = Tmat_avg(1);
T2_Matlab = Tmat_avg(2);

% --- Set outputs --------------------------------------------------------

trnOutputs(1) = T1_Matlab;
trnOutputs(2) = T2_Matlab;

% Tell TRNSYS we reached the end of the m-file without errors
mFileErrorCode = 0; 
return  % Do not update outputs at this call
