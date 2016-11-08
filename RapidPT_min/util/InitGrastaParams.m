function [ options, opts, opts2, status ] = InitGrastaParams( maxRank, iter, V )
%InitGrastaParams Initialize the parameters for GRASTA
%   Initialize the parameters for GRASTA. GRASTA is the library used for
%   matrix completion. Note: opts is used later..

    options.RANK         =      maxRank; %% (from inputs)
    options.rho          =      2;   
    options.ITER_MAX     =      iter; %% (from inputs)  
    options.ITER_MIN     =      options.ITER_MAX;   
    options.USE_MEX      =      0;    
    options.TOL          =      1e-8; 
    options.DIM_M        =      V; %% (from inputs)    

    opts                 =      struct(); 

    opts2.RHO            =      options.rho;  
    opts2.TOL            =      options.TOL;
    opts2.MAX_ITER       =      options.ITER_MAX;       

    status               =      struct(); 
    status.init          =      0; 
    
    fprintf('Grasta Parameters: ');
    disp(options);
    
end

