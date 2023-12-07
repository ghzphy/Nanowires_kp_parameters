% Input : params  : [mt,ml,k0,E0]
%         kpoints : [:,3]
%         Nx, Ny
% Output: eneriges: [:,30] Only return the smallest 30 eigenvalues.
function energies = E_K_reduce_001(params,kpoints)
    global N_bands Lx Ly
    [k_num,~] = size(kpoints);
    a = 3.80993; % hbar^2/(2*m0)
    mt_inv = a/params(1); 
    ml_inv = a/params(2); 
    M_inv = mt_inv - a;
    k0 = params(3);
    U = params(4) + ml_inv*k0^2; 
   
    KL = 100; dz = 1; %Ly = 30; Lz = 30;
    Ham00 = zeros(2*KL,2*KL);
    Ham01 = zeros(2*KL,2*KL);
    
    Nx = 100; Ny = 100;
    dis = zeros(Nx*Ny,1);
    num = zeros(Nx*Ny,2);
    num_new = zeros(KL,2);
    for i = 1:Nx
        for j = 1:Ny
            index = (i-1)*Ny + j;
            dis(index) = ((i-1)^2 + (j-1)^2)^0.5;
            num(index,:) = [i,j];
        end
    end
    [~,index] = sort(dis);
    for i =1: KL
        num_new(i,:) = num(index(i),:);
    end
    
    for k = 1: KL
        p = num_new(k,1);
        q = num_new(k,2);
        kp = p*pi/Lx;kq = q*pi/Ly;
        for k1 = 1: KL
            p1 = num_new(k1,1);
            q1 = num_new(k1,2);
            kp1 = p1*pi/Lx;kq1 = q1*pi/Ly;
            
            % H11 = 
            H00_11 = (ml_inv*2/(dz^2) + mt_inv*(kp^2 + kq^2) + U)*delta(p,p1)*delta(q,q1);
            % H22 = 
            H00_22 = (ml_inv*2/(dz^2) + mt_inv*(kp^2 + kq^2) + U)*delta(p,p1)*delta(q,q1);
            
            % H12 = -2*kx*ky/M
            if (q == q1) || (p == p1)
                H00_12 = 0;
            else
                H00_12 = 2*M_inv*(4*kp1/pi*p/(p^2-p1^2)*delta_odd(p,p1))*(4*kq1/pi*q/(q^2-q1^2)*delta_odd(q,q1));
            end
            H00_21 = H00_12;
            Ham00(2*k-1:2*k,2*k1-1:2*k1)= [H00_11,H00_12;H00_21,H00_22];
            
            % H11 = 
            H01_11 = (-ml_inv/(dz^2) + k0*1i*ml_inv/dz)*delta(p,p1)*delta(q,q1);
            % H22 = 
            H01_22 = (-ml_inv/(dz^2) - k0*1i*ml_inv/dz)*delta(p,p1)*delta(q,q1);
            % H12 = -2*kx*ky/M
            H01_12 = 0;
            H01_21 = H01_12;
            Ham01(2*k-1:2*k,2*k1-1:2*k1)= [H01_11,H01_12;H01_21,H01_22];
        end
    end
    
    Ham10 = Ham01';
    energies = zeros(k_num,N_bands);
    for i = 1:k_num
        kz = kpoints(i,3);
        Ham = Ham00 + exp(1j*kz)*Ham01 + exp(-1j*kz)*Ham10;
        vals =  eig(Ham);
        tmp = sort(real(vals));
        energies(i,:) = tmp(1:N_bands);
    end
end