% bounds is bounds_from_Lmask2bounds_function
function nuclei=Lbounds2nuclei(bounds)
nuclei=[];
for i=1:length(bounds)
    nuclei{i}=[bounds(i).r' bounds(i).c'];
end
end