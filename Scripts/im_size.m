function [IM1] = im_size(IM,k)
sz = size(IM);
for i = 1:sz(1)*k
    for j = 1:sz(2)*k
        IM1(i,j,:) = IM(ceil(i/k),ceil(j/k),:);
    end
end
end

