function [ string ] = stringForDATfile(plotType, bubble, SUPG, groupedFE )

if strcmp(plotType,'potential')
    string='dat_files/potL2Error';
elseif strcmp(plotType,'electronConc')
    string='dat_files/elecL2Error';
elseif strcmp(plotType,'holeConc')
    string='dat_files/holeL2Error';
end

if groupedFE
    string=strcat(string, 'trap');
else
    string=strcat(string, 'GQ');
end

if SUPG
    string=strcat(string, 'SUPG');
elseif bubble
    string=strcat(string, 'bubble');
end

string=strcat(string, '.dat');

end

