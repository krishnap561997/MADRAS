function output = encode_array(input,dir1_max, dir2_max, dir3_max, opts)

switch opts
    case 1   %Encode Data
        input = reshape(input, [numel(input)./3 , 3]);
        i1 = input(:,1);
        i2 = input(:,2);
        i3 = input(:,3);
        output = (i1-1) + dir1_max.*(i2-1) + dir1_max.*dir2_max.*(i3 - 1);
    case 2   %Decode Data
        input = reshape(input, [numel(input) , 1]);
        i3 = floor(input./dir1_max./dir2_max);
        rem1 = input - i3.*dir1_max.*dir2_max;
        i2 = floor(rem1./dir1_max);
        i1 = rem1 - i2.*dir1_max;    
        output = [(i1+1)';(i2+1)';(i3+1)']';
end



