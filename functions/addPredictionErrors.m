%% Prediction errors

if predErrStd > 0
    boudRatio = 1;
    if ~exist('pd','var')
        pd = 'normal';
    end
    
    % interactive workload
    a_err = random(pd, 0, mean(a)*predErrStd, size(a)) ;
   
    a = a_real - a_err; % a_pred
    a = max(a, 0); 
    a = min(a, mean(a_real)+boudRatio*mean(a_real)); 
        
    % batch jobs
    BS_err = random(pd, 0, mean(BS)*predErrStd, size(BS)) ;

    BS = abs(BS_real - BS_err); % BS_pred
    BS = max(BS, 0); 
    BS = min(BS, mean(BS_real)+boudRatio*mean(BS_real)); 
    
    % renewable 
    R_err  = random(pd, 0, mean(R)*predErrStd, size(R)) ;
    
    R = abs(R_real - R_err); % R_pred
    R = max(R, 0); 
    
    % prices
    GP_err  = random(pd, 0, mean(GP)*predErrStd, size(GP)) ;
    GP = abs(GP_real - GP_err); % GP_pred
    GP = max(GP, 0); 
    GP = min(GP, mean(GP_real)+boudRatio*mean(GP_real));
    
    % natural prices
    CRC_temp = CRC(2:T+1);
    CRC_err  = random(pd, 0, mean(CRC_temp)*predErrStd, size(CRC_temp)) ;
    CRC_temp = abs(CRC_real(2:T+1) - CRC_err); % GP_pred
    CRC_temp = max(CRC_temp, 0); 
    CRC_temp = min(CRC_temp, mean(CRC_real(2:T+1))+boudRatio*mean(CRC_real(2:T+1)));
    CRC(2:T+1) = CRC_temp;
    
    for y=1:N_y
        GP_Array(:,y)  = GP_Array(:,y) - random(pd, 0, mean(GP_Array(:,y))*predErrStd, size(GP_Array(:,y)));        
        GP_Array(:,y) = max(GP_Array(:,y), 0); 
        GP_Array(:,y) = min(GP_Array(:,y), mean(GP_Array_real(:,y))+boudRatio*mean(GP_Array_real(:,y)));
        
        a_Array(:,y)  = a_Array(:,y)  - random(pd, 0, mean(a_Array(:, y))*predErrStd, size(a_Array(:, y)));
        a_Array(:,y)  = max(a_Array(:,y), 0); 
        a_Array(:,y)  = min(a_Array(:,y), mean(a_Array_real(:,y))+boudRatio*mean(a_Array_real(:,y)));
        
        BS_Array(:,y)  = BS_Array(:,y) - random(pd, 0, mean(BS_Array(:, y))*predErrStd, size(BS_Array(:, y)));
        BS_Array(:,y)  = max(BS_Array(:,y), 0); 
        BS_Array(:,y)  = min(BS_Array(:,y), mean(BS_Array_real(:,y))+boudRatio*mean(BS_Array_real(:,y)));
        
        CRC_Array(2:T+1,y) = CRC_Array(2:T+1,y) - random(pd, 0, mean(CRC_Array(2:T+1,y))*predErrStd, size(CRC_Array(2:T+1,y)));        
        CRC_Array(2:T+1,y) = max(CRC_Array(2:T+1,y), 0); 
        CRC_Array(2:T+1,y) = min(CRC_Array(2:T+1,y), mean(CRC_Array_real(2:T+1,y))+boudRatio*mean(CRC_Array_real(2:T+1,y)));        
    end
else
    N_samples = 1;
end
