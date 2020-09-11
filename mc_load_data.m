function data = mc_load_data(Template,subs)

N = length(subs);

Subject = subs{1};
tmp = mc_load_datafile(mc_GenPath(Template));
s = size(tmp);
s = [N s];
data = zeros(s);
data(1,:,:) = tmp;

for i = 2:N
    Subject = subs{i};
    tmp = mc_load_datafile(mc_GenPath(Template));
    data(i,:,:) = tmp;
end
