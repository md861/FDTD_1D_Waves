%Created for the purpose of research and fun!
%Features:
%1D - scalar wave equation
%Backward Euler in time and Central Difference in space
%Periodic boundary condition
%NB: This is not an optimised code
%Author --- M.Drolia 
%Written -- July-2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define Spatial parameters
l_s = 100;                %Length of space
nel_s = 1000;             %# of elements in space
nnod_s = nel_s + 1;       %# of nodes in space
bcnodes_s = {1,nnod_s};   %indices of spatial boundaries
del_x = l_s/nel_s;        %spatial step-size

%Define Temporal parameters
l_t = 100;              %Length of time
nel_t = 200;            %# of elements in time
nnod_t = nel_t + 1;     %# of nodes in time
icnodes_t = {1,2};      %indices of time boundaries (initial condition)
del_t = l_t/nel_t;      %temporal step-size

%Define problem parameters
c = 1;                  %phase velocity of the wave equation


%Generate space-time mesh, solution vector, rhs vector and system matrix
x = [0:del_x:l_s];              %space variable
t = [0:del_t:l_t];              %time variable
f = zeros(nnod_t,nnod_s);       %Variable that stores field values
sys_mat = zeros(nnod_s,nnod_s); %System matrix
rhs_vec = zeros(nnod_s,1);      %RHS vector


%Build system matrix and source vector   
  %Update initial conditions
  n = icnodes_t{1};          %temporal index
  for k = 1:nnod_s
      f(n,k) = ic_term(t(n),k,nnod_s);
  end
  n = icnodes_t{2};          
  for k = 1:nnod_s
      f(n,k) = ic_term(t(n),k,nnod_s);
  end
  %Implement d2f/dt2 - c2d2f/dx2 = source_term
  c2_g = (del_x*del_x)/(del_t*del_t); %Grid speed square - for brevity
  term_offd = -c*c/c2_g;              %Off-diagonal constants
  term_diag = 1-(2*term_offd);        %Diagonal constants
  %at bc1 - Periodic boundary condition
  k = bcnodes_s{1};
  sys_mat(k,bcnodes_s{2}-1) = term_offd;
  sys_mat(k,k) = term_diag;
  sys_mat(k,k+1) = term_offd;
  %at bc2 - Periodic boundary condition
  k = bcnodes_s{2};
  sys_mat(k,k-1) = term_offd;
  sys_mat(k,k) = term_diag;
  sys_mat(k,bcnodes_s{1}+1) = term_offd;
  %populate rest of the terms in the system matrix
  for k = 2:nnod_s-1
      sys_mat(k,k-1) = term_offd;
      sys_mat(k,k) = term_diag;
      sys_mat(k,k+1) = term_offd;
  end 
  


%Solve for the scalar field for each time
for n = icnodes_t{2}+1:nnod_t
    %Build the RHS vector = src + previous time solutions    
    for k = 1:nnod_s
        rhs_vec(k) = (del_t*del_t*src_term(t(n),k,nnod_s)) + ...
                     2*f(n-1,k) - f(n-2,k);
    end
    f(n,:) = sys_mat\rhs_vec;
end

%Plot
h_fig = figure(1);
for n = 1:nnod_t
    plot(x,f(n,:))
    xlim([0 l_s]);
    ylim([0 0.3])
    title('Sinusoidal response for the wave equation');
    xlabel('Space: 1-Dimensional');
    ylabel('Scalar field');
    M(n) = getframe;
end
%movie(M,1)


%Function definitions 
    %Source term
    function src_val = src_term(t,x_index,nnod_s)
        %Sinusoidal source -----
        l_temp = 10;   
        lx_temp = 5;
        if (x_index == ceil(nnod_s/lx_temp))
            src_val = sin(3*pi*t/l_temp);
        else
            src_val = 0;
        end
        %Source free ---------
        %src_val = 0;
    end
    %Initial condition
    function ic_val = ic_term(t,x_index,nnod_s)
        %Sinusoidal
        ic_val = src_term(t,x_index,nnod_s);
        %Delta function ------
        %lx_temp = 5;
        %if (x_index == ceil(nnod_s/lx_temp))
        %    ic_val = 1;
        %else
        %    ic_val = 0;
        %end
    end
    