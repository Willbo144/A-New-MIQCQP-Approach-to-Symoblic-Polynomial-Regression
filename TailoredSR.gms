Sets
i 'data point i'
p 'predictor variable p e.g. x(i, p)'
                    
n /1*7/

xt(n) /1*3/

nterm(n)

layers /1*2/
xlayers(layers,xt) /1. 2*3
                    2. 1
                    /

* Operators including binary, uniray and predictors/constants
b /'+','x1','x2','*','const','cube'/
* List of avaliable variables that you dont want to use
nb(b)//
*List of binary variables
bi(b) /'+','*'/
*list of uniary variables
uni(b) /'cube'/
*constants and predictors
con(b) /'x1','x2','const'/
*predictors only
data(b) /'x1','x2'/
*Symmetric operators
sym(b) /'+','*'/

alias(i,j);

nterm(n) = not xt(n);

Variables
v(i,n) * value of nodes
vl(i,xt) * value of a nodes left child
vr(i,xt) *value of a nodes right child
er * error
cst(n) * value of constant 
;
Binary Variables
y(b,n) * binary variable for all operators.
; 
Parameters
z(i)
std
vup2(i,n)  * final upper bound for each node
vlo2(i,n)  * final lower bound for each node
Mup(i,n,b) * big M up
Mlo(i,n,b) * big M lo
Msc(i,xt)  * big M for symmetry cuts
vup(i,n,b) * intermediate upper limit for the value of a node
vlo(i,n,b) * intermediate lower limit for the value of a node
x(i,p) data points for two predictors x1 and x2
* Upper limit for node value
Vlim /1.0e4/
* Node limit
Nlim /7/
*Safety Margin
s /1/
*Upper and lower constant limits
conUP /10/
conLO 
;
* imports predictor data from csv
$call csv2gdx predictor_data.csv id=x index=1 values=2..lastCol useheader=y
$gdxIn predictor_data.gdx
$load i = Dim1
$load p = Dim2
$load x
$gdxIn
conLO = -conUP;
* Function 
z(i) =  2*x(i,'1')/exp(0.1*x(i,'2'));
* Standards deviation of response data
std = sqrt(sum(i,sqr(z(i)-sum(j,z(j))/card(z)))/(card(z)-1));

*sets terminal node limits
vup2(i,nterm) = s+max(smax(p, x(i, p)),conUP);
vlo2(i,nterm) = -s+min(smin(p, x(i, p)),conLO);

* sets node limits for layers above terminal nodes
loop((xlayers(layers,xt(n))),

 vup(i,n,'+') =  s+vup2(i,n+ord(n))+vup2(i,n+(ord(n)+1));
 vlo(i,n,'+') = -s+vlo2(i,n+ord(n))+vlo2(i,n+(ord(n)+1));


 vup(i,n,'*') = s+max(vup2(i,n+ord(n))*vup2(i,n+(ord(n)+1)),vup2(i,n+ord(n))*vlo2(i,n+(ord(n)+1)),vlo2(i,n+ord(n))*vup2(i,n+(ord(n)+1)),vlo2(i,n+ord(n))*vlo2(i,n+(ord(n)+1)));
 vlo(i,n,'*') = -s+min(vup2(i,n+ord(n))*vup2(i,n+(ord(n)+1)),vup2(i,n+ord(n))*vlo2(i,n+(ord(n)+1)),vlo2(i,n+ord(n))*vup2(i,n+(ord(n)+1)),vlo2(i,n+ord(n))*vlo2(i,n+(ord(n)+1)));


 vup(i,n,'x1') = s+smax(j,x(j,'1'));
 vlo(i,n,'x1') =  -s+smin(j,x(j,'1'));
 vup(i,n,'x2') = s+smax(j,x(j,'2'));
 vlo(i,n,'x2') =  -s+smin(j,x(j,'2'));

 vup(i,n,'const') = s+ conUP;
 vlo(i,n,'const') = -s+ conLO;
 
 vup(i,n,'cube') =  s +power(vup2(i,n+ord(n)),3);
 vlo(i,n,'cube') = -s + power(vlo2(i,n+ord(n)),3);

 
vup2(i,n) = min(Vlim,smax(b,vup(i,n,b)));
vlo2(i,n) = max(-Vlim,smin(b,vlo(i,n,b)));


);

vup2(i,'1') = z(i)+2*std;
vlo2(i,'1') = z(i)-2*std;

* calculates big M values
Mup(i,n,b) = vup(i,n,b) + abs(vlo2(i,n));
Mlo(i,n,b) = vlo(i,n,b) - vup2(i,n);
Msc('1',xt(n)) = vlo2('1',n+ord(n))-vup2('1',n+(ord(n)+1));

display Mup,Mlo,Msc,vup2,vlo2;

Equations
error,
plusUP,plusLO,
timesUP,timesLO,
cubeUP,cubeLO,
constUP,constLO,constUP2,constLo2,constantUP,
op,
terminal,terminal2,terminal3,terminal4,terminal5,terminal6,terminal7,terminal8,
node0up,node0lo,nodetup,nodetlo,
child1,child2,
binar,binar2,
pred,pred2,predmin,
num,Termin,
active,
uniray,uniray2,
nodelimit,
operatorlim,
symmetrycut
;



* square error function
error..sum(i,sqr(z(i)-v(i,'1'))) =e= er;

* Equations for non terminal nodes
plusUP(i,xt)..vl(i,xt)+vr(i,xt)-v(i,xt)=l= Mup(i,xt,'+')*(1-y('+',xt));
plusLO(i,xt)..vl(i,xt)+vr(i,xt)-v(i,xt)=g= Mlo(i,xt,'+')*(1-y('+',xt));
terminal(i,xt).. x(i,'1')-v(i,xt)=l= Mup(i,xt,'x1')*(1-y('x1',xt));
terminal2(i,xt).. x(i,'1')-v(i,xt)=g= Mlo(i,xt,'x1')*(1-y('x1',xt));
terminal3(i,xt).. x(i,'2')-v(i,xt)=l= Mup(i,xt,'x2')*(1-y('x2',xt));
terminal4(i,xt).. x(i,'2')-v(i,xt)=g= Mlo(i,xt,'x2')*(1-y('x2',xt));
timesUP(i,xt).. (vl(i,xt)*vr(i,xt))-v(i,xt) =l= Mup(i,xt,'*')*(1-y('*',xt));
timesLO(i,xt).. (vl(i,xt)*vr(i,xt))-v(i,xt) =g= Mlo(i,xt,'*')*(1-y('*',xt));
constUP(i,xt).. cst(xt)-v(i,xt) =l= Mup(i,xt,'const')*(1-y('const',xt));
constLo(i,xt).. cst(xt)-v(i,xt) =g= Mlo(i,xt,'const')*(1-y('const',xt));
cubeUP(i,xt).. power(vl(i,xt),3) - v(i,xt) =l= Mup(i,xt,'cube')*(1-y('cube',xt));
cubeLO(i,xt).. power(vl(i,xt),3) - v(i,xt) =g= Mlo(i,xt,'cube')*(1-y('cube',xt));



*Equations for terminal nodes
terminal5(i,nterm).. x(i,'1')-v(i,nterm)=l= Mup(i,nterm,'x1')*(1-y('x1',nterm));
terminal6(i,nterm).. x(i,'1')-v(i,nterm)=g= Mlo(i,nterm,'x1')*(1-y('x1',nterm));
terminal7(i,nterm).. x(i,'2')-v(i,nterm)=l= Mup(i,nterm,'x2')*(1-y('x2',nterm));
terminal8(i,nterm).. x(i,'2')-v(i,nterm)=g= Mlo(i,nterm,'x2')*(1-y('x2',nterm));
constUP2(i,nterm).. cst(nterm)-v(i,nterm) =l= Mup(i,nterm,'const')*(1-y('const',nterm));
constLo2(i,nterm).. cst(nterm)-v(i,nterm) =g= Mlo(i,nterm,'const')*(1-y('const',nterm));

* defines limits for non-terminal node values
node0up(i,xt).. v(i,xt) =l= vup2(i,xt)*sum(b,y(b,xt));
node0lo(i,xt).. v(i,xt) =g= vlo2(i,xt)*sum(b,y(b,xt));
*defines limits for terminal node values
nodetup(i,nterm).. v(i,nterm) =l= vup2(i,nterm)*sum(con,y(con,nterm));
nodetlo(i,nterm).. v(i,nterm) =g= vlo2(i,nterm)*sum(con,y(con,nterm));

* sets limit for the constant between 10 and -10 to prevent it being very very small
constantUP(n).. abs(cst(n)) =l= conUP*y('const',n);


* only one operator per node
op(n).. sum(b,y(b,n))=l=1;

* specifies parent child relations 
child1(i,xt(n)).. vl(i,xt) =e= v(i,n+ord(n));
child2(i,xt(n)).. vr(i,xt) =e= v(i,n+(ord(n)+1));
* if a binary operator is used then both children must be active 
binar(xt(n)).. sum(bi,y(bi,xt)) =l= sum(b,y(b,n+ord(n)));
binar2(xt(n)).. sum(bi,y(bi,xt)) =l= sum(b,y(b,n+(ord(n)+1)));

* if a uniray operator is used then only the left hand node can be active and its cant be a constant i.e. 2^2 or Log 3 etc
uniray(xt(n)).. sum(uni,y(uni,xt)) =l= sum(b,y(b,n+ord(n)))-y('const',n+ord(n));
uniray2(xt(n))..sum(uni,y(uni,xt)) =l= 1-sum(b,y(b,n+(ord(n)+1)));

* if a parent is a constant or a predictor then its children must be inactive
pred(xt(n)).. sum(con,y(con,xt))=l=1-sum(b,y(b,n+ord(n)));
pred2(xt(n))..sum(con,y(con,xt))=l=1-sum(b,y(b,n+(ord(n)+1)));

* if a parent is inactive then its children must be inactive
active(xt(n)).. sum(b,y(b,n+ord(n))) + sum(b,y(b,n+(ord(n)+1))) =l= 2*sum(b,y(b,xt));

* no pair of children can be constants only or x1 x1 , x2 x2
num(xt(n)).. y('const',n+ord(n))+ y('const',n+(ord(n)+1)) =l= 1;

*Ensures when a symmetric operator is used the left hand node is larget than the right
symmetrycut(xt).. vl('1',xt)-vr('1',xt) =g= Msc('1',xt)*(1-sum(sym,y(sym,xt)));

* Solution must use at least one predictor 
predmin.. sum((data,n),y(data,n)) =g= 1;
*Terminal nodes cannot assume binary or uniary operators
Termin.. sum((bi,nterm),y(bi,nterm))+ sum((uni,nterm),y(uni,nterm)) =e= 0;
*Node limiter
nodelimit.. sum((b,n),y(b,n)) =l= Nlim;
*Sets binary variables for undesired operators to 0
operatorlim.. sum((nb,n),y(nb,n)) =e= 0;
   

Model equ /all/;
Option minlp=BARON;
Solve equ using minlp minimizing er;
display x, z, cst.l, y.l;



