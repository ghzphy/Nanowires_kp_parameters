% Input : params  : [L,M,N,E0]
%         kpoints : [:,3]
% Output: eneriges: [:,3]
function energies = bulk_valence(params,kpoints)
    [k_num,~] = size(kpoints);
    a = 3.80993; % hbar^2/(2*m0)
    L = a*params(1); 
    M = a*params(2);
    N = a*params(3);
    E0 = params(4); 
    
    energies = zeros(k_num,3);
    for i = 1:k_num
        H = zeros(3);
        kx = kpoints(i,1);
        ky = kpoints(i,2);
        kz = kpoints(i,3);
        H(1,1) = L*kx^2 + M*(ky^2 + kz^2) + E0;
        H(2,2) = L*ky^2 + M*(kx^2 + kz^2) + E0;
        H(3,3) = L*kz^2 + M*(kx^2 + ky^2) + E0;
        H(1,2) = N*kx*ky;H(1,3) = N*kx*kz;
        H(2,1) = N*kx*ky;H(2,3) = N*ky*kz;
        H(3,1) = N*kx*kz; H(3,2) = N*ky*kz;
        energies(i,:) = eig(H);
    end
end