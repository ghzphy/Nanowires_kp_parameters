function kpath = kpoints2kpath(kpoints)
   [k_num,~] = size(kpoints);
   kpath = zeros(1,k_num);
   tmp = 0;
   for i = 2:k_num
       tmp = tmp + ((kpoints(i,1)-kpoints(i-1,1))^2 + (kpoints(i,2)-kpoints(i-1,2))^2 + (kpoints(i,3)-kpoints(i-1,3))^2)^0.5;
       kpath(i) = tmp;
   end
end
