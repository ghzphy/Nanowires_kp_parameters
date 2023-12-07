% Input : params  : [L,M,N,E0]
%         kpoints : [:,3]
% Output: eneriges: [:,N_E] Only return the N_E eigenvalues.
function energies = H_K_reduce_110(params,kpoints)
    global N_bands Lx Ly
    [k_num,~] = size(kpoints);
    a = 3.80993; % hbar^2/(2*m0)
    L = a*params(1); M = a*params(2);N = a*params(3); E0 = params(4);
    
    KL = 100; dz = 1; %Lx = 30; Ly = 30;
    Ham00 = zeros(3*KL,3*KL);
    Ham01 = zeros(3*KL,3*KL);
    
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
            
            H00_11 = ((L+M)/(dz^2) + M*kp^2 + (L+M)/2*kq^2 + E0)*delta(p,p1)*delta(q,q1);
            H00_22 = ((L+M)/(dz^2) + L*kp^2 + (L+M)/2*kq^2 + E0)*delta(p,p1)*delta(q,q1);
            H00_33 = (M*2/(dz^2) + L*kp^2 + M*kq^2 + E0)*delta(p,p1)*delta(q,q1);
            % H12=N*kx*ky; 
            H00_12 = N/2*(2/(dz^2) - kq^2)*delta(p,p1)*delta(q,q1); 
            H00_21 = H00_12;
            
            % H13=N*kx*kz;%H23=N*ky*kz
            if (q == q1) || (p == p1)
                H00_13 = 0;
                H00_23 = 0;
            else
                H00_13 = N/sqrt(2)*(4*kp1/pi*p/(p^2-p1^2)*delta_odd(p,p1))*(4*kq1/pi*q/(q^2-q1^2)*delta_odd(q,q1));
                H00_23 = -N/sqrt(2)*(4*kp1/pi*p/(p^2-p1^2)*delta_odd(p,p1))*(4*kq1/pi*q/(q^2-q1^2)*delta_odd(q,q1)); 
            end
            H00_31 = H00_13;
            
            H00_32 = H00_23;
            Ham00(3*k-2:3*k,3*k1-2:3*k1)= [H00_11,H00_12,H00_13;H00_21,H00_22,H00_23;H00_31,H00_32,H00_33];
            
            if q == q1
                H01_11 = -(L+M)/(2*dz^2)*delta(p,p1)*delta(q,q1);
                H01_22 = -(L+M)/(2*dz^2)*delta(p,p1)*delta(q,q1);
            else
                H01_11 = -(L+M)/(2*dz^2)*delta(p,p1)*delta(q,q1) + (L-M)/(2*dz)*(4*kq1/pi*q/(q^2-q1^2)*delta_odd(q,q1))*delta(p,p1);
                H01_22 = -(L+M)/(2*dz^2)*delta(p,p1)*delta(q,q1) - (L-M)/(2*dz)*(4*kq1/pi*q/(q^2-q1^2)*delta_odd(q,q1))*delta(p,p1);
            end

            H01_33 = -M/(dz^2)*delta(p,p1)*delta(q,q1);
            %H12=N*kx*ky; H13=N*kx*kz;
            H01_12 = -N/(2*dz^2)*delta(p,p1)*delta(q,q1);
            H01_21 = H01_12;
            
            %H23=N*ky*kz
            if p ~= p1
                H01_13 = -N/(sqrt(2)*2*dz)*(4*kp1/pi*p/(p^2-p1^2)*delta_odd(p,p1))*delta(q,q1);
                H01_23 = -N/(sqrt(2)*2*dz)*(4*kp1/pi*p/(p^2-p1^2)*delta_odd(p,p1))*delta(q,q1);
            else
                H01_13 = 0;
                H01_23 = 0;
            end
            H01_31 = H01_13;
            H01_32 = H01_23;
            Ham01(3*k-2:3*k,3*k1-2:3*k1)= [H01_11,H01_12,H01_13;H01_21,H01_22,H01_23;H01_31,H01_32,H01_33];
        end
    end
    
    Ham10 = Ham01';
    energies = zeros(k_num,N_bands);
    for i = 1:k_num
        kz = kpoints(i,3);
        Ham = Ham00 + exp(1j*kz)*Ham01 + exp(-1j*kz)*Ham10;
        vals =  eig(Ham);
        tmp = sort(real(vals),'descend');
        energies(i,:) = tmp(1:N_bands);
    end
end