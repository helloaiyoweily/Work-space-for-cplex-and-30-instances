function data_reading(F,Cdata)
Cdata_axis=Cdata(:,2:3);
Cdata_demands=Cdata(:,4:5);
D=F(find(F(:,end)==1),:);       
P=F(find(F(:,end)==1),:);       
F=F;       
I=Cdata(find(Cdata(:,end)==1),:);     
J=Cdata(find(Cdata(:,end)==2),:);     
S=Cdata(find(Cdata(:,end)==3),:);     
U=Cdata;      
V=10;       
K=5;    
cusNum=size(Cdata,1);
DNum=size(D,1);
PNum=size(P,1);
SNum=size(S,1);
cusI=size(I,1);
cusJ=size(J,1);
cusS=size(S,1);
vNum=V;
kNum=K;
facNum=size(F,1);
Qh=Cdata(:,5);     
Gh=Cdata(:,6)  ;     
eh1=Cdata(:,6)  ;     
eh2=Cdata(:,7)  ;     
ef1=F(:,6)  ;      
ef2=F(:,7)  ;      
drh =[]  ;       
dwf =[]  ;       
trh =[]  ;        
e=20  ;      
l=30  ;    
T=52  ;     
Tv=10  ;      
Tk=15  ;      
fv=0.06  ;      
fk=0.1  ;      
Pv=6.4  ;      
Pk=6.4  ;      
Cf=1500  ;      
Cv=200  ;      
Ck=600  ;      
Mv=20000  ;      
Mk=25000 ;      
Gf=F(:,5)  ;     
CCf=F(:,6)  ;     
If=F(:,7)  ;      
tmof=[]  ;      
tmfo=[]  ;      
Tmo=[]  ;      
ndvrh=[]  ;    
npvrh=[]  ;      
dtvf=[]  ;     
atvh=[]  ;      
LB=[]  ;      
UB=[]  ;      
RT=7200  ;     
GAP_AL=0  ;     
GAP_UL=0  ;        
save ParametersAll 
end





