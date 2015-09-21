function cleanupWaitbar
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);