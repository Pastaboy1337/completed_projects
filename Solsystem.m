clear; clc;
%%ger ett slumpmässigt system om detta önskas

RandomNumberPlanet = randi([3,10],1);
RandomPlanetRadius = randi([25,120],1,RandomNumberPlanet);
RandomPlanetTime = randi([20,30],1,RandomNumberPlanet);
RandomNumberMoons = randi([3,9],1);
RandomPlanetMoonAttachment = randi([1,RandomNumberPlanet],1,RandomNumberMoons);
RandomMoonOrbitTime = randi([5,60],1,RandomNumberMoons);
RandomMoonRadius = randi([2,15],1,RandomNumberMoons);
R = RandomPlanetRadius;
T = RandomPlanetTime;
M1 =RandomPlanetMoonAttachment;
M2 =RandomMoonOrbitTime;
M3= RandomMoonRadius;
deltaT = 0.1;

%% Huvudprogrammet, tar in värden och skickar dem till main funktionen,
%sparar sedan ner den inspelade filmen som en .avi fil på datorn i hela
% 30 fps.stänger sedan alla öppna fönster så att dem öppnas igen nästa gång
% programmet körs.
film = main(R,T,M1,M2,M3,deltaT);
v = VideoWriter("film");
v.open;
v.writeVideo(film);
v.close;
movie(film,1,30)
close all
%% funktionen där allt sker, tar in värden på följande sätt:
% R = radien för planeten i obeatämd längdenhet
% T = omloppstiden för planeten i obestämd tidsenhet,
% simulationen följer 2 varv av den långsammaste planeten
% M1 = vilken av planeterna i arrayen som månen kretsar kring.
% M2 = planeternas omloppsbana obestämd tidsenhet
% M3 = Planeternas omkretsradie i obestämd längdenhet
% deltaT = tidssteget på varje itteration av simulationen,
% lägre värde = mer noggrant resultat, kräver dock mer ram-minne samt tid.
% OBS! kör inte stora värden på T med små värden på deltaT, programmet
% kommmer eventuellt fylla ramminnet
function film = main(R,T,M1,M2,M3,deltaT)
arguments
    % ser till att endast numeriska värden kommer in, förhindrar att programmet
    % crashar pga felaktiga inputs.
    R,T,M1,M2,M3,deltaT (1,:) {mustBeNumeric, mustBePositive}
end
%% kollar så att alla månar sitter på en planet samt att alla planet och
% mån-arraysen har samma storlek så att allt är definierat. .
if max(M1)>length(R)||length(R)~=length(T)||length(M1)~=length(M2)...
        ||length(M2)  ~= length(M3) ||length(M1) ~= length(M3)

    switch true
        case max(M1)>length(R)
            [~,pos] = max(M1);
            ErrorMessage = "Moon must be attached to a planet, value at index " + int2str(pos)+...
                " in M1 must be smaller than/or equal to " + int2str(length(R));
        case length(R)~=length(T)
               ErrorMessage = "All planetary arrays must have the same size"
        case length(M1)~=length(M2)
                ErrorMessage = "All Moon arrays must have the same size"
        case length(M2)  ~= length(M3)
                   ErrorMessage = "All planetary arrays must have the same size"
        case length(M1) ~= length(M3)
            ErrorMessage = "All planetary arrays must have the same size"
    end
    msgbox(ErrorMessage);
    return
end
%% sätter/ hämtar värden som ska användas i kalkyler
DeltaT = deltaT;
PlanetRadies = R;
PlanetOrbitTimes = T;
numPlanets = length(PlanetRadies);
planetCords = zeros(1,1,1);
MaxT = 2*max(PlanetOrbitTimes);
Maxr = max(R)+max(M3);
MoonRadies = M3;
MoonOrbitTimes= M2;
PlanetConnection = M1;
numMoons = length(MoonRadies);
FrameNumber = 0;
PlanetColour = rand([numPlanets, 3]);
MoonColour = rand([numMoons, 3]);
RandomStarNumber = randi([100,150],1);
RandomStarSize = randi([1,9],1,RandomStarNumber);
k = round(Maxr *1.1,0);
RandomStarCoordinatesX = randi([-k, k],1,RandomStarNumber );
RandomStarCoordinatesY = randi([-k, k],1,RandomStarNumber );
time = -DeltaT;
%% skapar två tomma anottations som visar hur lång tid det har gått på
% simulationer samt hur lång tid det ör kvar av simulatinen, värdena på
% dessa sätts senare i den stora loopen.
msg ='';
msg2='';
timeCounter = annotation("textbox",[.85 .6 .1 .2],'String',msg,'EdgeColor','none');
timeCounterRemaining = annotation("textbox",[.85 .4 .1 .2],'String',msg2,'EdgeColor','none');
%% beräknar alla planeternas omloppsbanor
for i = 1:1:numPlanets
    angleIndex = 0;
    for alfa = linspace(0,2*pi,1E4)
        angleIndex = angleIndex+1;
        x = cos(alfa)*PlanetRadies(i);
        y = sin(alfa)*PlanetRadies(i);
        planetCords(1,angleIndex,i) = x;
        planetCords(2,angleIndex,i) = y;
    end
end
%% beräknar vinkelhastigheten för alla planeter
AngularVelocity = [];
for i = 1:1:numPlanets
    AngularVelocity(i) = (2*pi)/PlanetOrbitTimes(i); %#ok<AGROW>
end
for i = 1:1:length(MoonOrbitTimes)
    MoonAngularVelocity(i) = (2*pi)/MoonOrbitTimes(i);
end
%% Skapar en array med ett slumpmässigt antal stjärnor av olika storlek på
% slumpmssiga platser. och ritar ut dem.
% ritar samtidigt ut alla planeternas omloppsbanor.
for i =1:1:RandomStarNumber
    ax = gca;
    set(gca,'Color','k','xtick',[],'ytick',[])
    ax.YLim = [-1.1.*Maxr 1.1.*Maxr];   %
    ax.XLim = [-1.1.*Maxr 1.1.*Maxr];
    axis equal
    hold on
    scatter(RandomStarCoordinatesX(i),RandomStarCoordinatesY(i),RandomStarSize(i),'*','w');
    scatter(0,0,250,'filled','y');
    for PlanetNmber = 1:1:numPlanets
        set(gca,'Color','k','xtick',[],'ytick',[],'PlotBoxAspectRatio',[1 1 1])
        plot(planetCords(1,:,PlanetNmber),planetCords(2,:,PlanetNmber),'-',color = PlanetColour(PlanetNmber,1:3))
    end
    %%
end
%% Sparar ner denna plot som en .png bild på datorn och läsar sedan in den
% under namnet I
f= gcf;
exportgraphics(f,'background.png','Resolution',600);
I = imread('background.png');

%% tar bort den tidigare plotens värden så att de inte fyller minnet mer än
% nödvändigt.
clear [RandomStarSize, RandomStarCoordinatesY, RandomStarCoordinatesX,...
    RandomStarNumber];
%% PlotLoopen, här sker själva simulationen.
for T = 0:DeltaT:MaxT
    background = image(xlim,ylim,I); % importerar bakgrunden
    uistack(background,'bottom'); % sätter den längs bak i plotten
    hold on
    PlanetPositions = []; % tömmer arrayen
    ax = gca;
    %fixar axlarna så att ploten inte visar mer än nödvändigt.
    %fixar även aspect ration så att plotten blir en kub
    ax.YLim = [-1.1.*Maxr 1.1.*Maxr];
    ax.XLim = [-1.1.*Maxr 1.1.*Maxr];
    axis equal
    set(gca,'Color','k','PlotBoxAspectRatio',[1 1 1])

    %% räknar ut planeternas position
    for PlanetIndex = 1:1:numPlanets
        X = cos(AngularVelocity(PlanetIndex)*T)*PlanetRadies(PlanetIndex);
        Y = sin(AngularVelocity(PlanetIndex)*T)*PlanetRadies(PlanetIndex);
        PlanetPositions(1,1,PlanetIndex) = X; %#ok<AGROW>
        PlanetPositions(1,2,PlanetIndex) = Y; %#ok<AGROW>
    end
    %% plottar ut alla planeternas position
    for PlanetNmber = 1:1:numPlanets
        scatter(PlanetPositions(:,1,PlanetNmber),PlanetPositions(:,2,PlanetNmber),75,'filled',MarkerFaceColor=PlanetColour(PlanetNmber,1:3))
    end
    %%
    MoonPositions = [];
    %%räknar ut månarnas position
    for MoonIndex = 1:1:numMoons
        X = cos(MoonAngularVelocity(MoonIndex)*T)*MoonRadies(MoonIndex) ;
        Y = sin(MoonAngularVelocity(MoonIndex)*T)*MoonRadies(MoonIndex);
        MoonPositions(1,1,MoonIndex) = X; %#ok<AGROW>
        MoonPositions(1,2,MoonIndex) = Y; %#ok<AGROW>
        MoonPositions(:,:,MoonIndex) = MoonPositions(:,:,MoonIndex) + PlanetPositions(:,:,PlanetConnection(MoonIndex)); %#ok<AGROW>
    end
    %% Räknar ut månardnas omloppsbanor
    for MoonNumber = 1:1:numMoons
        MoonCordsComplete =[1,1,1];
        MoonCords = [1,1,1];

        angleIndex = 0;
        %%
        for alfa = linspace(0,2*pi,360)
            angleIndex = angleIndex+1;
            % ger x och y kordinater i cirkel
            x = cos(alfa)*MoonRadies(MoonNumber);
            y = sin(alfa)*MoonRadies(MoonNumber);
            MoonCords(1,angleIndex,MoonNumber) = x;
            MoonCords(2,angleIndex,MoonNumber) = y;
        end
        %% räknar ut kring vilken punkt månarnas omkrets har sin mittpunkt
        for j = 1:1:size(MoonCords,2)
            MoonCordsComplete(1,j,MoonNumber)= MoonCords(1,j,MoonNumber) + PlanetPositions(:,1,PlanetConnection(MoonNumber));
            MoonCordsComplete(2,j,MoonNumber)= MoonCords(2,j,MoonNumber) + PlanetPositions(:,2,PlanetConnection(MoonNumber));
        end
        %%
        %% ritar ut månarnas position samt deras kretsbanor.
        plot(MoonCordsComplete(1,:,MoonNumber),MoonCordsComplete(2,:,MoonNumber),'-',color = MoonColour(MoonNumber,1:3));
        scatter(MoonPositions(:,1,MoonNumber),MoonPositions(:,2,MoonNumber),20,'filled',MarkerFaceColor=MoonColour(MoonNumber,1:3))

    end
    %% fixar annotationen så att den displaryar rätt värde
    time = time+DeltaT;
    timeremaining = MaxT - time;
    msg = sprintf('time passed: %f days', time);
    set(timeCounter,'string',msg);
    msg2 = sprintf('time remaining: %f days',timeremaining);
    set(timeCounterRemaining,'string',msg2);
    FrameNumber = FrameNumber+1;
    film(FrameNumber) = getframe; %#ok<AGROW>
    hold off;
end
%%
end
