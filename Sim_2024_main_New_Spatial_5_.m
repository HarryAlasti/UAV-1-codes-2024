%% Assumptions:
% 1. The spatial signal is distributed over a coordinated, known area.
% 2. The spatial signal distribution varies slowly.
% 3. The UAVs have suitable sensors, and GPS for localization purposes.
% 4. There is an accessible information fusioning center (FC) for the
% required processing.

%% Purpose:
% Monitoring the spatial signal over time based on data and communication
% efficient approach.

clc; close all; clear all; 
pause(1);

%% Simulation assumptions:
n_std = 0.0;                                                        % Standard deviation of the additive white noise to sensor's observation.
Nlevels0 = 3;                                                       % Initial assumption for the number of levels.
Npnt = 8;
step = 1;
Grid_pnt = 10; 
Error = [];
Error_real = [];
dist = [];
jmp0 = 1;
Delta = [];
Level_flag = 0;                                                     % Level_flag = 1 for equally-spaced levels, otherwise optimally-spaced.
SG_flag = 1;                                                         % SG_flag = 1, SG is applied,     SG_flag = 0, SG not applied
%% Generation of the spatial signal
name = sprintf('NewWave_%d',1);
load(name);
v0 = 0 : step : 100;
[X , Y] = meshgrid(v0);
Z = SpatialSignal(X,Y,1, n_std);
%figure(1000); mesh(X,Y,Z);
%pause(1);
Lmin = min(Z(:));
Lmax = max(Z(:));
sp0 = Lmax - Lmin;

%% Getting a pre-estimation of the spatial signal using UAV observations
% Based on the two rough estimations (that can be combined into one), the
% UAV system capture the signal strength at a number of points, that are
% assumed over the grid.

X1 = linspace(0,100,Grid_pnt);
Y1 = X1;
X1 = [X1 X1];
Y1 = [Y1 Y1(end:-1:1)];
Z1 = SpatialSignal(X1(:),Y1(:),1,n_std);
dist_0 = 2 * sqrt(2) * 100;

Zold = griddata(X1(:), Y1(:), Z1(:), X, Y, 'v4');
%figure(10); mesh(X,Y,Zold);
%pause;
%pause;

Nlevels = Nlevels0;
Lmax = max(Zold(:));
Lmin = min(Zold(:));

if (Level_flag == 1)
    Levels = linspace(Lmin, Lmax, Nlevels+2);
    Levels = Levels(2:end-1);
else
    [p,xx] = hist(Zold(:), 500);
    p = p / sum(p);
    Levels = Lloyd_Max_2(p, xx, Nlevels)';
end

delta = min(abs(Levels(1:end-1) - Levels(2:end)));
Delta = [Delta, delta];

LEVELS = Levels;
AccLevels = LEVELS;

Space = [];
JMP = [];
SP = [];
AccDist = 0;
NLevels = [Nlevels0];
spr = 1;

for npnt = 1 : Npnt
    Cx = contourc(v0, v0, Zold, LEVELS);
    Cx0 = contourc(v0, v0, Z, LEVELS);
    [cx , Dist, Frstpnt, Num_pnts] = Parser3(Cx);
    [cx0 , Dist0, Frstpnt0, Num_pnts0] = Parser3(Cx0);
    [sz1 , sz2] = size(Frstpnt);

    for count_cnt = 1 : sz1
        pnt = Frstpnt(count_cnt , :);        
        [mini , num] = min(sqrt((cx0(:,1) - pnt(1)).^2 + (cx0(:,2) - pnt(2)).^2 + (cx0(:,3) - pnt(3)).^2));
        cont_pce = cx0(num, 4);
        jj = find(cx0(: , 4) == cont_pce);
        x = cx0(jj , 1:3);
        Space = [Space ; x];
        AccDist = AccDist + Dist0(cont_pce);
    end    
    dist = [dist , AccDist];
    
    Zx = griddata(Space(:,1), Space(:,2), Space(:,3), X,Y,'v4');

    Error = [Error, mean(abs(Zx(:) - Zold(:)))];
    Error_real = [Error_real, mean(abs(Zx(:) - Z(:)))];

    Zold = Zx;

   if (SG_flag == 1 && npnt > 1)
        %jmp =  ceil(jmp *(1 - 2*abs(Error(npnt-1) - Error(npnt))/(Error(npnt) + Error(npnt-1)))); 

        jmp =  jmp + ceil((1 + 2*abs(Error(npnt-1) - Error(npnt))/(Error(npnt) + Error(npnt-1)))); 
        delta = delta * abs(1 - 2*(Error(npnt-1) - Error(npnt))/(Error(npnt) + Error(npnt-1)));
   elseif(SG_flag == 1 && npnt == 1)
       jmp = 1;
   end

  if (SG_flag == 0)
        jmp = jmp0;
   end

   Delta = [Delta , delta];
    JMP = [JMP , jmp];

    if (npnt < Npnt)
        Nlevels = Nlevels0 + jmp;
        Nlevels0 = Nlevels;
        Lmax = max(Zx(:));
        Lmin = min(Zx(:));
       
        if (Level_flag == 1)
            Levels = linspace(Lmin, Lmax, Nlevels+2);
            Levels = Levels(2:end-1);
        else
            [p,xx] = hist(Zx(:),200);
            p = p /sum(p);     
            Levels = Lloyd_Max_2(p, xx, Nlevels);
        end
    end   
    sp = Lmax - Lmin;
    spr = sp/sp0;
    SP = [SP , spr];
%{
  Levels'
  AccLevels
  delta
  pause
%}
   
    if (SG_flag == 1 && npnt > 1)    
        ThisLevels = [];
        for lvl = 1 : length(Levels)
            [Min,ind] = min(abs(AccLevels - Levels(lvl)));
            if (Min >= delta)
               ThisLevels = [ThisLevels ,  Levels(lvl)];
            end
        end    
    else
        ThisLevels = Levels';
    end

   Nend = NLevels(end);
   NLevels = [NLevels , length(ThisLevels) + Nend];
   
    if (npnt < Npnt)         
        AccLevels = [AccLevels , (ThisLevels(:))'];

        if (npnt == 1)
            LEVELS = Levels';
        else
            LEVELS = ThisLevels;
        end 
        %Nlevels = length(ThisLevels);
    end

end

dist = dist_0 + dist;
%figure(2000);  plot(NLevels(1:Npnt), (Error(1:Npnt)), 'd-r',NLevels(1:Npnt), (Error_real(1:Npnt)), 's-g'); xlabel('# of Levels'); ylabel('Mean Absolute Error '); legend("Learning Error", "Actual Error");grid on;
figure(2500);  plot(NLevels(1:Npnt), 20*log10(Error(1:Npnt)), 'r',NLevels(1:Npnt), 20*log10(Error_real(1:Npnt)), 'g'); xlabel('# of Levels'); ylabel('Mean Absolute Error (dB)'); legend("Learning Error + SG", "Actual Error+SG");grid on;
figure(3000);  plot(NLevels(1:Npnt), dist(1:Npnt), 'd-b'); xlabel('# of Levels'); ylabel('UAV paced distance '); grid on;
figure(2700);  plot(NLevels(1:Npnt), SP(1:Npnt), 'd-b'); xlabel('# of Levels'); ylabel('Signal span '); grid on;
figure (2900); plot(dist(1:Npnt), Error_real(1:Npnt), 'd-b'); xlabel('Cost (Flying Distance)'); ylabel('MAE (dB)'); grid on;

if (SG_flag == 1)
    figure(3500); plot(1:Npnt,Delta(1:Npnt)); grid on;
end

%NLevels
%pause

figure(3100); mesh(X,Y,abs(Zx - Z)); grid on;


%% 
clc; clear all; 
pause(1);

%% Simulation assumptions:
n_std = 0.0;                                                        % Standard deviation of the additive white noise to sensor's observation.
Nlevels0 = 3;                                                       % Initial assumption for the number of levels.
Npnt = 13;
step = 4;
Grid_pnt = 25; 
Error = [];
Error_real = [];
dist = [];
jmp0 = 1;
Delta = [];
Level_flag = 0;                                                     % Level_flag = 1 for equally-spaced levels, otherwise optimally-spaced.
SG_flag = 0;                                                         % SG_flag = 1, SG is applied,     SG_flag = 0, SG not applied
%% Generation of the spatial signal
name = sprintf('NewWave_%d',1);
load(name);
v0 = 0 : step : 100;
[X , Y] = meshgrid(v0);
Z = SpatialSignal(X,Y,1, n_std);
%figure(1000); mesh(X,Y,Z);
Lmin = min(Z(:));
Lmax = max(Z(:));
sp0 = Lmax - Lmin;

%% Getting a pre-estimation of the spatial signal using UAV observations
% Based on the two rough estimations (that can be combined into one), the
% UAV system capture the signal strength at a number of points, that are
% assumed over the grid.

X1 = linspace(0,100,Grid_pnt);
Y1 = X1;
X1 = [X1 X1];
Y1 = [Y1 Y1(end:-1:1)];
Z1 = SpatialSignal(X1(:),Y1(:),1,n_std);
dist_0 = 2 * sqrt(2) * 100;

Zold = griddata(X1(:), Y1(:), Z1(:), X, Y, 'v4');
%figure(10); mesh(X,Y,Zold);

Nlevels = Nlevels0;
Lmax = max(Zold(:));
Lmin = min(Zold(:));

if (Level_flag == 1)
    Levels = linspace(Lmin, Lmax, Nlevels+2);
    Levels = Levels(2:end-1);
else
    [p,xx] = hist(Zold(:), 500);
    p = p / sum(p);
    Levels = Lloyd_Max_2(p, xx, Nlevels)';
end

delta = min(abs(Levels(1:end-1) - Levels(2:end)))/2;
Delta = [Delta, delta];

LEVELS = Levels;
AccLevels = LEVELS;

Space = [];
JMP = [];
SP = [];
AccDist = 0;
NLevels = [Nlevels0];

for npnt = 1 : Npnt
    Cx = contourc(v0, v0, Zold, LEVELS);
    Cx0 = contourc(v0, v0, Z, LEVELS);
    [cx , Dist, Frstpnt, Num_pnts] = Parser3(Cx);
    [cx0 , Dist0, Frstpnt0, Num_pnts0] = Parser3(Cx0);
    [sz1 , sz2] = size(Frstpnt);

    for count_cnt = 1 : sz1
        pnt = Frstpnt(count_cnt , :);        
        [mini , num] = min(sqrt((cx0(:,1) - pnt(1)).^2 + (cx0(:,2) - pnt(2)).^2 + (cx0(:,3) - pnt(3)).^2));
        cont_pce = cx0(num, 4);
        jj = find(cx0(: , 4) == cont_pce);
        x = cx0(jj , 1:3);
        Space = [Space ; x];
        AccDist = AccDist + Dist0(cont_pce);
    end    
    dist = [dist , AccDist];
    
    Zx = griddata(Space(:,1), Space(:,2), Space(:,3), X,Y,'v4');

    Error = [Error, mean(abs(Zx(:) - Zold(:)))];
    Error_real = [Error_real, mean(abs(Zx(:) - Z(:)))];

    Zold = Zx;

   if (SG_flag == 1 && npnt > 1)
        jmp =  ceil(jmp *(1 + 2*abs(Error(npnt-1) - Error(npnt))/(Error(npnt) + Error(npnt-1)))); 
        delta = delta * abs(1 - 2*(Error(npnt-1) - Error(npnt))/(Error(npnt) + Error(npnt-1)));
   elseif(SG_flag == 1 && npnt == 1)
       jmp = 1;
   end

  if (SG_flag == 0)
        jmp = jmp0;
   end

   Delta = [Delta , delta];
    JMP = [JMP , jmp];

    if (npnt < Npnt)
        Nlevels = Nlevels0 + jmp;
        Nlevels0 = Nlevels;
        Lmax = max(Zx(:));
        Lmin = min(Zx(:));
       
        if (Level_flag == 1)
            Levels = linspace(Lmin, Lmax, Nlevels+2);
            Levels = Levels(2:end-1);
        else
            [p,xx] = hist(Zx(:),200);
            p = p /sum(p);     
            Levels = Lloyd_Max_2(p, xx, Nlevels);
        end
    end   
    sp = Lmax - Lmin;
    spr = sp/sp0;
    SP = [SP , spr];
   
    if (SG_flag == 1 && npnt > 1)    
        ThisLevels = [];
        for lvl = 1 : length(Levels)
            [Min,ind] = min(abs(AccLevels - Levels(lvl)));
            if (Min >= delta)
               ThisLevels = [ThisLevels ,  Levels(lvl)];
            end
        end   
    end
    if ((SG_flag == 1 && npnt == 1) || (SG_flag == 0))
        ThisLevels = Levels';
    end

   Nend = NLevels(end);
   NLevels = [NLevels , length(ThisLevels) + Nend];

    if (npnt < Npnt)         
        AccLevels = [AccLevels , (ThisLevels(:))'];

        if (npnt == 1)
            LEVELS = Levels';
        else
            LEVELS = ThisLevels;
        end 
        %Nlevels = length(ThisLevels);
    end
end
dist = dist_0 + dist;

%figure(2000); hold on; plot(NLevels(1:Npnt), (Error(1:Npnt)), 'd--r',NLevels(1:Npnt), (Error_real(1:Npnt)), 's--g'); xlabel('# of Levels'); ylabel('Mean Absolute Error '); legend("Learning Error", "Actual Error");grid on;
figure(2500);  hold on;  plot(NLevels(1:Npnt), 20*log10(Error(1:Npnt)), '--r',NLevels(1:Npnt), 20*log10(Error_real(1:Npnt)), '--g'); xlabel('# of Levels'); ylabel('Mean Absolute Error (dB)'); legend("Learning Error + SG", "Actual Error + SG", "Learning Error - No SG", "Actual Error - No SG");grid on;
figure(3000);  hold on; plot(NLevels(1:Npnt), dist(1:Npnt), 'd--b'); xlabel('# of Levels'); ylabel('UAV paced distance '); grid on;
figure(2700);  hold on; plot(NLevels(1:Npnt), SP(1:Npnt), 'd--b'); xlabel('# of Levels'); ylabel('Signal span '); grid on;
figure (2900); hold on; plot( dist(1:Npnt), Error_real(1:Npnt),'s--r'); xlabel('Cost (Flying Distance)'); ylabel('MAE (dB)'); grid on;
if (SG_flag == 1)
    figure(3500); hold on; plot(1:Npnt,Delta(1:Npnt),'--r'); grid on;
end
figure(3200); mesh(X,Y,abs(Zx - Z)); grid on;