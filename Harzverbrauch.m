% Laminatdicke given in datasheet is with resin!
% Not exactly the thickness of the textile!

layerArea = 1; % m2
layerThickness = 0.308*10^-3; % m
layerVolume = layerArea*layerThickness;


layerMass = 500; % g/m2
fvc = 0.35;


resinMass = 220; % g/m2
fiberMass = layerMass - resinMass;


fiberVolume = fvc*layerVolume;
resinVolume = layerVolume - fiberVolume;

fiberThickness = fiberVolume/layerArea;

%%
fvc2 = 0.4;
resinDensity2 = 1.151*10^6; % g/m3

fiberVolume2 = fiberVolume;
layerVolume2 = fiberVolume2/fvc2;
resinVolume2 = layerVolume2 - fiberVolume2;
resinMass2 = resinDensity2*resinVolume2;

Harzverbrauch_LW_fvc40 = resinMass2;


%%

