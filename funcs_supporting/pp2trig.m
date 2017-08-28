

function trig = pp2trig (pp,thresh)

if nargin < 2
    thresh = 1000;
end

ppd = diff(pp);
ind = [1, find(ppd > thresh)+1];
ind = ind(ind < length(pp));        % Get rid of anything that's too long
trig = pp(ind);

end