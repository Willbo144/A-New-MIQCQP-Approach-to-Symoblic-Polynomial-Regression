Sets
i
p
b /con,x1,x2,x1negone,x2negone,x1two,x2two,x1negtwo,x2negtwo,logx1,logx2,lnx1,lnx2,x1three,x2three,x1negthree,x2negthree
   x1half,x2half,x1onehalf,x2onehalf,x1negonehalf,x2negonehalf,x1neghalf,x2neghalf,x1twohalf,x2twohalf,x1negtwohalf,x2negtwohalf/;


Parameters
z(i)
x(i,p);
Binary variables
y(b);
Variables
C(b)
er
f
fit(i);

$call csv2gdx predictor_data.csv id=x index=1 values=2..lastCol useheader=y
$gdxIn predictor_data.gdx
$load i = Dim1
$load p = Dim2
$load x
$gdxIn

z(i) = power(x(i,'1'),3)/(6.5+x(i,'1'));
C.l(b) = 0;
Equations
error,obj,limits1,fitted,termlim;
termlim.. sum(b,y(b)) =l= 20;
fitted(i).. C('con')+
           (C('x1half')*vcpower(x(i,'1'),0.5))+(C('x2half')*vcpower((x(i,'2')),0.5))+
           (C('x1')*x(i,'1'))+(C('x2')*(x(i,'2')))+
           (C('x1onehalf')*vcpower(x(i,'1'),1.5))+(C('x2onehalf')*vcpower(x(i,'2'),1.5))+
           (C('x1two')*sqr(x(i,'1')))+(C('x2two')*sqr(x(i,'2')))+
           (C('x1twohalf')*vcpower(x(i,'1'),2.5))+(C('x2twohalf')*vcpower(x(i,'2'),2.5))+
           (C('x1three')*power(x(i,'1'),3))+(C('x2three')*power(x(i,'2'),3))+
           (C('x1neghalf')*vcpower(x(i,'1'),-0.5))+(C('x2neghalf')*vcpower((x(i,'2')),-0.5))+
           (C('x1negone')*(1/x(i,'1')))+(C('x2negone')*(1/(x(i,'2'))))+
           (C('x1negonehalf')*vcpower(x(i,'1'),-1.5))+(C('x2negonehalf')*vcpower(x(i,'2'),-1.5))+
           (C('x1negtwo')*power(x(i,'1'),-2))+(C('x2negtwo')*power(x(i,'2'),-2))+
           (C('x1negtwohalf')*vcpower(x(i,'1'),-2.5))+(C('x2negtwohalf')*vcpower(x(i,'2'),-2.5))+
           (C('x1negthree')*power(x(i,'1'),-3))+(C('x2negthree')*power(x(i,'2'),-3))+
           (C('lnx1')*log(x(i,'1')))+(C('lnx2')*log(x(i,'2')))+
           (C('logx1')*log10(x(i,'1')))+(C('logx2')*log10(x(i,'2'))) =e= fit(i);
*THIS IS THE VALUE OF THE ERROR           
error.. sum(i,sqr(z(i)-fit(i)))=e= f;
*THIS IS NOT THE VALUE OF THE ERROR 
obj.. (card(z)*(log10(sum(i,sqr(z(i)-fit(i)))/card(z))))+(2*sum(b,y(b))*log10(log10(card(z)+1)))=e= er;
limits1(b).. abs(C(b)) =l= 100*y(b); 
Model equ /all/;
Option minlp=BARON;
Solve equ using minlp minimizing er;
display x,z,fit.l,C.l;
option decimals =8; display f.l;
