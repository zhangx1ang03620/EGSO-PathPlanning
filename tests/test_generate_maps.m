clc; clear; close all;
[Maps, Map_Info] = GenerateEGSOMaps();
print(gcf, fullfile(tempdir, 'EGSO_Test_Maps'), '-dpng', '-r200');
close(gcf);
fprintf('\nSaved %s\n', fullfile(tempdir, 'EGSO_Test_Maps.png'));
