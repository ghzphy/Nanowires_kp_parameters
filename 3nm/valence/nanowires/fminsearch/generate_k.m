function kpoints = generate_k(begin,end_k,len)
    kpoints = zeros(len,3);
    for i = 1:len
        kpoints(i,1) = begin(1)+(i-1)*(end_k(1)-begin(1))/(len-1);
        kpoints(i,2) = begin(2)+(i-1)*(end_k(2)-begin(2))/(len-1);
        kpoints(i,3) = begin(3)+(i-1)*(end_k(3)-begin(3))/(len-1);
    end
end