function [sorted]=sortpop (unsorted,column)
[T,I]=sort(unsorted(:,column));
sorted=unsorted(I,:);