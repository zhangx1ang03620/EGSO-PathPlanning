% Application of simple limits/bounds
function s = Bounds( s, Xmin, Xmax)

% Keep the position within bounds
index = find(s > Xmax);
s(index) = randi(Xmax);
index = find(s < Xmin);
s(index) = randi(Xmax);
s = round(s);
