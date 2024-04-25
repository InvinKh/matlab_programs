function [nd_arr] = derivativeKV(a_arr, l)
nd_arr=zeros(numel(a_arr(:,1)), numel(a_arr(1,:)));
switch l
    case 1
        d_arr = diff(a_arr, 1, 2);
        for i=1:numel(a_arr(:,1))
            for j=1:numel(a_arr(1,:))-1
                nd_arr(i,j)=nd_arr(i,j) + d_arr(i,j)/2;
                nd_arr(i,j+1) = nd_arr(i,j+1) + d_arr(i,j)/2;
            end
        end
    case 2
        d_arr = diff(a_arr, 1, 1);
        for j=1:numel(a_arr(1,:))
            for i=1:numel(a_arr(:,1))-1
                nd_arr(i,j)=nd_arr(i,j) + d_arr(i,j)/2;
                nd_arr(i+1,j) = nd_arr(i+1,j) + d_arr(i,j)/2;
            end
        end
end