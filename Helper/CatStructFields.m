function S = CatStructFields(dim, varargin)

F = cellfun(@fieldnames,varargin,'uni',0);
assert(isequal(F{:}),'All structures must have the same field names.')
T = [varargin{:}];
S = struct();
F = F{1};

for k = 1:numel(F)
    S.(F{k}) = cat(dim,T.(F{k}));
end

new_sess = false(size(S.mouseID));
for i=1:numel(T)
   if i==1
       new_sess(1) = true;
       counter=1+length(T(1).mouseID);
   else
       new_sess(counter) = true;
       counter=counter+length(T(i).mouseID);
   end
end

S.new_sess = new_sess;