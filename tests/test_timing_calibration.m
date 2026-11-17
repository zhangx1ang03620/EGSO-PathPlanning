clc; clear;
[Maps, Map_Info] = GenerateEGSOMaps();
close(gcf);

G = Maps{5};  % 90x90, largest/densest map
Xmax = size(G,1); X_min = 1; dim = Xmax - 2;
fobj = @(x) fitness(x, G);

names = {'EGSO', 'ARO', 'INFO', 'SBOA', 'GWO', 'ACO'};
for k = 1:numel(names)
    fn = str2func(names{k});
    tic;
    fn(50, 100, X_min, Xmax, dim, fobj, G);
    elapsed = toc;
    fprintf('%-8s (N=50,T=100, map 90x90): %.2f s\n', names{k}, elapsed);
end
