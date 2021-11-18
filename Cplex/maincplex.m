clear all;  close all; clc; tic
tstart=tic;
load F;
load Cdata;
data_parameter_reading(F,Cdata); load parametersAll
keyvariables(Cdata,F); load variables1;
options = cplexoptimset('cplex');
options.timelimit=RT;
D=pdist(Cdata_axis);
Drh=squareform(D);
trh=Drh/4;
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
TC=0;TCV=0;TCK=0;TCF=0
for w=1
    for o=1
        for i=1:vNum
            
            for r=1:facNum
                for h=1:cusNum
                    TC=TC+ xvrh*fv*Pv*Drh(r,h);
                end
            end
            
            for r=1:cusNum
                for h=1:cusNum
                    TC=TC+ xvrh*fv*Pv*Drh(r,h);
                end
            end
            
            for r=1:cusNum
                for h=1:facNum
                    TC=TC+ xvrh*fv*Pv*Drh(r,h);
                end
            end
        end
        
        for v=1:vNum
            for r=1:facNum
                for h=1:cusNum
                    TCV =TCV+xvrh(v,r,h)*Mv/T
                end
            end
        end
        
        for v=1:vNum
            for r=1:facNum
                for h=1:cusNum
                    TCV=TCV+xvrh(v,r,h)*(e*max((eh1(h)-atvh(v,h)),0) + l*max((atvh(v,h)-eh2(h)),0))
                end
            end
        end
    end
    for p=1
        for k=1:kNum
            for o=1:facNum
                for w=1:facNum
                    for f=1:facNum
                        TCK=TCK+fk*Pk*xkowf(k,o,w,f)*Drh(w,f)
                    end
                end
            end
          end
        for o=1:facNum
            TCK=TCK+Mk*max(Tmo(o)/Ck)
        end       
    end 
    for i=1:facNum
         TCF=TCF+F(i,6)+F(i,7)-zf(i)*F(i,8);
    end
end
TC=TCV+TCK+TCF;
f=TC;
FC=[];
for h=1:cusNum
    for v=1:vNum
        for f=1:facNum
            FC=[FC;sum(xvhf(v,v,f))==1]; 
        end
    end
end
for v=vNum
    for r=1:facNum
        for h=1:cusNum
           FC=[FC;sum(xvrh(v,:,:))==1];
        end
    end
end
    for v=1:vNum
       for r=1:cusNum
           FC=[FC;sum(xvrh(v,r,:))==sum(xvrh(v,:,r))];      
       end 
    end
for r=1:facNum
   for n=1:facNum
       for v=1:vNum
           FC=[FC;sum(xvrh(v,r,:))==sum(xvrh(v,:,r))];                   
       end       
   end
end
for o=1:facNum
    for k=1:kNum
        FC=[FC; sum(xkowf(k,o,o,:))==sum(xkowf(k,o,:,o))];
    end
end
for o=1:facNum
    for k=1:kNum
        for o=1:facNum
            for f=1:facNum
                FC=[FC; sum(xkowf(k,o,:,f))==sum(xkowf(k,o,f,:))];
            end
            FC=[FC; sum(xkowf(k,o,o,:))==sum(xkowf(k,o,:,o))];
        end
        
    end
end
for k=1:kNum
    for o=1:facNum
        for f=1:facNum
            FC=[FC; sum(xkowf(k,o,:,f))==sum(xkowf(k,o,f,:))];
        end
    end
end
for k=1:kNum
    for o=1:facNum
        for f=1:facNum
            FC=[FC; sum(xkowf(k,o,:,f))==sum(xkowf(k,o,f,:))];
        end
    end
end
for v=1:vNum
    for r=1:cusNum
        for h=1:cusI
            FC=[FC;sum(ndvrh(v,r,h)+npvrh(v,r,h)+xvrh(v,r,h)*(Gh(h)-Qh(h)))<=Cv];
        end
        for hh=(cusI+cusS+1):cusNum
            FC=[FC;sum(ndvrh(v,r,hh)+npvrh(v,r,hh)+xvrh(v,r,hh)*(Gh(hh)-Qh(hh)))<=Cv];
        end
    end
end
for k=1:kNum
    for o=1:DNum
        FC=[FC;sum(xkof(k,o,:).*tmof(o,:))<=Ck];
    end
end
for k=1:kNum
    for o=(DNum+1):facNum
        FC=[FC;sum(xkof(k,o,:).*tmfo(:,o))<=Ck];
    end
end
for o=1:DNum
    for f=1:DNum
        for h=1:DNum
            tmof(o,f)=tmof(o,f)+yhof(h,o,f)*Qh(h);
        end
        for h=1:SNum
            tmof(o,f)=tmof(o,f)+yhof(h,o,f)*Qh(h);
        end
    end
end
for o=1:PNum
    for f=1:PNum
        for h=1:PNum
            tmfo(f,o)=tmfo(f,o)+yhof(h,f,o)*Qh(h);
        end
        for h=SNum
            tmfo(f,o)=tmfo(f,o)+yhof(h,f,o)*Gh(h);
        end
    end
end
for w=1
    for o=1:DNum
        Tmo(o)=sum(tmof(o,:)) ;
    end
    for o=(DNum+1):facNum
        Tmo(o)=sum(tmfo(:,o));
    end
end
for f=1:DNum
    for h=1:(cusI)
   FC=[FC;(sum(xvhf(:,h,f)*Qh(h))+sum(yhof(h,:,f)*Qh(h)))<=Cf(f)]; 
    end
     for h=(cusI+cusJ+1):cusNum
   FC=[FC;(sum(xvhf(:,h,f)*Qh(h))+sum(yhof(h,:,f)*Qh(h)))<=Cf(f)]; 
      end
end
for f=(DNum+1):facNum
     FC=[FC;( sum(xvhf(v,(cusI+1):cusNum,f)*Gh(h))+sum(yhof((cusI+1):cusNum,(DNum+1):facNum,f)*(Gh(h)) ))<=Cf];
end
for v=1:vNum
    for f=1:facNum
        FC=[FC;dtvf(v,f)<=ef1(f)];
        FC=[FC;dtvf(v,f)>=ef2(f)];
    end
end
for v=1:vNum
    for r=1:cusNum
        for h=1:cusNum
            FC=[FC;(((atvh(v,r)+trh(r,h)))*xvrh(v,r,h))>=(eh1(h)*xvrh(v,r,h))];
            FC=[FC;(((atvh(v,r)+trh(r,h)))*xvrh(v,r,h))<=(eh2(h)*xvrh(v,r,h))];        
        end
    end
end
for v=1:vNum
    for r=1:cusNum
        FC=[FC;sum(xvrh(v,r,:).*trh(r,:))<=Tv]
    end
end
for v=1:vNum
    for f=1:facNum
        FC=[FC;sum(xvrh(v,:,:))<=(sum(xvhf(v,:,f))-1)];
    end
end
for k=1:kNum
    for o=1:DNum
        FC=[FC;(sum(xkowf(k,o,:,:)))<=(sum(xkof(k,o,:))-1)];
    end
end
for k=1:kNum 
    for o=1:PNum
        FC=[FC;(sum(xkowf(k,o,:,:)))<=(sum(xkof(k,o,:))-1)];
    end
end
ops = sdpsettings( 'solver','cplex');
sol = solvesdp(FC,f,ops);
f   = double(f);
xvhf  = double(xvhf);
xvrh  = double(xvrh);
yhof  = double(yhof);
xkowf = double(xkowf);
zf    = double(zf);
tmof  = double( tmof );      	
tmfo  = double( tmfo );      	
Tmo   = double( Tmo );      	
ndvrh = double( ndvrh );      
npvrh = double( npvrh );       
dtvf  = double( dtvf );       
atvh  = double( atvh );     

 
