function keyvariables(Cdata,F)
data_parameter_reading(F,Cdata); load parametersAll
tmof  = sdpvar(facNum,facNum,'full');      
tmfo  = sdpvar(facNum,facNum,'full');      	
Tmo   = sdpvar(facNum,facNum,'full');      
ndvrh = sdpvar(vNum,cusNum,cusNum,'full');      
npvrh = sdpvar(vNum,cusNum,cusNum,'full');      
dtvf  = sdpvar(vNum,facNum,'full');       
atvh  = sdpvar(vNum,cusNum,'full');       
xvhf  = binvar(vNum,cusNum,facNum,'full') ;     
xvrh  = binvar(vNum,cusNum,cusNum,'full') ;      
yhof  = binvar(cusNum,facNum,facNum,'full') ;   
xkowf = binvar(kNum,facNum,facNum,facNum,'full') ;      
xkof  = binvar(kNum,facNum,facNum,'full') ;  
zf    = binvar(facNum,facNum,'full') ; 
save variables1
end