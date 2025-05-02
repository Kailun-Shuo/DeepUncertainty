%%%% Construction performance simulation with Simio (simulation platform can be replaced), Deep Uncertainty Analysis with PRIM

n=1;
ReplicationNo=1000;    
LHS=lhsdesign(ReplicationNo,9); % Latin hypercube sampling for intervals of deep uncertainty factors
LHS1=lhsdesign(ReplicationNo,9); % Latin hypercube sampling for Scale factors

for ii=1:ReplicationNo
    % Sampling scale factors with ranges from 1 to 5
    WP_WandCFormworkTime1_Scale(ii)=ceil(LHS1(ii,1)*5);
    WP_BandSFormworkTime1_Scale(ii)=ceil(LHS1(ii,2)*5);
    WP_CraneWorkingTime1_Scale(ii)=ceil(LHS1(ii,3)*5);
    WP_BandSRebarInstallTime1_Scale(ii)=ceil(LHS1(ii,4)*5);
    WP_ConcretePumpTime1_Scale(ii)=ceil(LHS1(ii,5)*5);
    PlatformRisingTime1_Scale(ii)=ceil(LHS1(ii,6)*5);
    WP_StiffcolumnHoistTime1_Scale(ii)=ceil(LHS1(ii,7)*5);
    WP_WRebarInstallTime1_Scale(ii)=ceil(LHS1(ii,8)*5);
    WP_CRebarInstallTime1_Scale(ii)=ceil(LHS1(ii,9)*5);
                             
    % Constructing new probability distributions from standard beta distribution
    % Generate random samples and select samples above the 75th percentile
    d0=betarnd(5,5,1,1e4);
    pd = makedist('Beta','a',5,'b',5);
    percentile_75 = icdf(pd,0.75);
    
    % WP_WandCFormworkTime1 distribution
    l1 = d0(d0>=percentile_75);
    l1 = repmat(l1, 1, WP_WandCFormworkTime1_Scale(ii));
    y1 = [d0 l1];  
    z11 = sort(y1(randperm(length(y1))));
    for i=1:ReplicationNo
        z111(i)=z11(round(length(z11)*(i/ReplicationNo)));
    end
    z111=[min(z11) z111]; % Constructed quantiles of the new distribution

    % Similarly for WP_BandSFormworkTime1
    d5=betarnd(5,5,1,1e4);
    l5 = d5(d5>=percentile_75);
    l5 = repmat(l5, 1, WP_BandSFormworkTime1_Scale(ii));
    y5 = [d5 l5];  
    z55 = sort(y5(randperm(length(y5))));
    for i=1:ReplicationNo
        z555(i)=z55(round(length(z55)*(i/ReplicationNo)));
    end
    z555=[min(z55) z555];

    % Similar steps are repeated for WP_CraneWorkingTime1, WP_BandSRebarInstallTime1, etc.
    % (code abbreviated here for brevity; same logic as above applies)
    
    % Generate WP_WandCFormworkTime based on selected intervals from LHS
    f=round(LHS(ii,1)*(ReplicationNo-1)+1);
    g1=z11(z111(f)<=z11 & z11<=z111(f+1));           
    WP_WandCFormworkTime{ii}=['Math.If(NoOfWandCFormworkWorker>30,',num2str((g1(unidrnd(length(g1)))*0.38+0.4)),',',num2str((g1(unidrnd(length(g1)))*0.43+0.3)),')']; 

    % WP_BandSFormworkTime similarly selected
    f=round(LHS(ii,2)*(ReplicationNo-1)+1);
    g5=z55(z555(f)<=z55 & z55<=z555(f+1));           
    WP_BandSFormworkTime{ii}=['Math.If(NoOfBandSFormworkWorker>30,',num2str((g5(unidrnd(length(g5)))*0.67+1.43)),',',num2str((g5(unidrnd(length(g5)))*0.94+1.33)),')']; 

    % Similar logic for other WP parameters (e.g., CraneWorkingTime, ConcretePumpTime, etc.)
    % (code abbreviated here for brevity)

    % Constants and predefined parameters for this simulation run
    NoOfWandCFormworkWorker='30';
    NoOfBandSFormworkWorker='30';
    CraneType='1';
    NoOfCrane='1';
    NoOfBandSRebarWorker='15';
    TimeofPath1='41';
    TimeofPath2='47';
    TimeofPath3='49';
    TimeofPath4='37';
    Concrete2='0';
    Concrete3='0';
    Concrete4='0';
    NoOfPump='2';
    NoOfRebarTruck='1';
    NoOfSteelColumnTruck='1';
    TruckType='1';
    RebarTransportTime='110';
    AlFrameworkRelated_WasteRate='0.04';
    Steel_WasteRate='0.077';
    Concrete_WasteRate='0.0899';
    NoOfLayer='10';
    NoOfCRebarWorker='30';
    NoOfWRebarWorker='30';
                 
    % Call simulation model (Simio or similar, the CallSimio is a API function, whose copyright is owned and can be downloaded at https://www.youtube.com/watch?v=_GpaW3MEgFo)
    [Time,TotalCost]=CallSimio(WP_WandCFormworkTime{ii},NoOfWandCFormworkWorker,WP_BandSFormworkTime{ii},NoOfBandSFormworkWorker,WP_CraneWorkingTime{ii},CraneType,NoOfCrane,NoOfBandSRebarWorker,WP_BandSRebarInstallTime{ii},TimeofPath1,TimeofPath2,TimeofPath3,TimeofPath4,Concrete2,Concrete3,Concrete4,NoOfPump,WP_ConcretePumpTime{ii},PlatformRisingTime{ii},WP_StiffcolumnHoistTime{ii},NoOfRebarTruck,NoOfSteelColumnTruck,TruckType,RebarTransportTime,AlFrameworkRelated_WasteRate,Steel_WasteRate,Concrete_WasteRate,NoOfLayer,WP_WRebarInstallTime{ii},WP_CRebarInstallTime{ii},NoOfCRebarWorker,NoOfWRebarWorker);
       
    % Record results
    Time_mean(ii)=nanmean(Time);
    Cost_mean(ii)=nanmean(TotalCost);
    
end

% Post-process results: calculate percentile values and deviations
T1(n,:)=Time_mean;
Time_mean1=sort(Time_mean);
Time_mean90(n,:)=Time_mean1(round(ReplicationNo*0.9));
Time_aver(n,:)=Time_mean1(round(ReplicationNo*0.5));
D1(n)=(Time_mean90(n,:)-Time_aver(n,:))/Time_aver(n,:);

C1(n,:)=Cost_mean;
Cost_mean1=sort(Cost_mean);
Cost_mean90(n,:)=Cost_mean1(round(ReplicationNo*0.9));
Cost_aver(n,:)=Cost_mean1(round(ReplicationNo*0.5));
D2(n)=(Cost_mean90(n,:)-Cost_aver(n,:))/Cost_aver(n,:);

% Selection of stable scenarios (adjust N accordingly)
N=1; % select scenario with better stability
X=[WP_WandCFormworkTime1_Scale' WP_BandSFormworkTime1_Scale' WP_CraneWorkingTime1_Scale' WP_BandSRebarInstallTime1_Scale' WP_ConcretePumpTime1_Scale' PlatformRisingTime1_Scale' WP_StiffcolumnHoistTime1_Scale' WP_WRebarInstallTime1_Scale' WP_CRebarInstallTime1_Scale'];     

% Binary outcome definition for PRIM analysis (Threshold set at 1710 for Time)
T11=T1(N,:);
Y1=T11'>=1710;

C11=C1(N,:);
Y2=C11'>=9075000;

% PRIM analysis parameters (PRIM toolbox's copyright is owned and can be download at http://www.cs.rtu.lv/jekabsons/regression.html#:~:text=Patient%20Rule%20Induction%20Method&text=PRIM%20is%20a%20method%20for,well%20as%20classification%2Dtype%20data.)
params = primparams(0.1,0.05,0.3);
[prim1, trajectory1] = primbuild(X, Y1, params,[0,0,0,0,0,0,0,0,0]);
[prim2, trajectory2] = primbuild(X, Y2, params,[0,0,0,0,0,0,0,0,0]);

% Print PRIM analysis results
primprint(prim1);
primprint(prim2);
