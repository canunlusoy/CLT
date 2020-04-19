function [fiberVolumeContent] = fvc(fiberVolume,matrixVolume)
fiberVolumeContent = fiberVolume/(fiberVolume+matrixVolume);
end

