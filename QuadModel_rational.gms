Sets
*Total number of terms
m /1*4/
*Data point d
d 
*Total numer of predictors
N(m) 
*Total number of quadratic terms
q(m)
alias(m,i)
alias(m,j);
Sets
* Mapping valid ranges of i and j onto the corrosponding value of q
Qi(m,i) / 2. 1*1
          3. 1*2
          4. 1*3/
Qj(m,j) /set.Qi/
* Double mapping to produce every possible valid combination of i and j for the corrosponding value of q
Qm(m,i,j) /2.1*1.1*1
           3.1*2.1*2
           4.1*3.1*3/
;

Parameters
x(d,N) *Input data
Y(d)  *Response data
B * Coefficient limit
BigM *Maximum number of terms;

* Imports predictor data from csv;
$call csv2gdx predictor_data.csv id=x index=1 values=1..lastCol useheader=y
$gdxIn predictor_data.gdx
$load d = Dim1
$load N = Dim2
$load x
$gdxIn

*Defines set of quadratic terms
q(m) = not N(m);
$offorder
Variables
a(m) Coefficients
be(m) 
er Error
a0 Inital coefficient
b0 
v(d,m) Value of terms
vb(d,m)
ans(d);
Binary Variables
z(m) *Binary variable controlling each term
zb(m)
wl(m,i) *Left component binary variable 
wr(m,j) * Right component binary variable 
gl(m,j)
gr(m,j);

B = 10000000;
BigM = 4;
*Sets values of predictor terms to the value of the corrosponding predictor at data point d
v.fx(d,N) = x(d,N);
vb.fx(d,N) = x(d,N);
*Initialises the program as a standard linear regression
z.l(q) = 0;
zb.l(m) = 0;
b0.fx = 1;
* Function to find
Y(d) = 1*exp(-(40700/8.314)*(1/(10*x(d,'1'))-1/(10*37.3)));
 

Equations
error,
Termlim,
Coefflim,
decompL,decompR,
values
Termlimb,
Coefflimb,
decompLb,decompRb,
valuesb
surrogate
;
* Numerator
error.. sum(d, sqr(Y(d) - ans(d))) =e= er;
Termlim.. sum(m, z(m)) =l= BigM;
Coefflim(m).. abs(a(m)) =l= z(m)*B;
decompL(q).. sum(Qi(q,i),wl(q,i)) =e= 1;
decompR(q).. sum(Qj(q,j),wr(q,j)) =e= 1; 
values(Qm(q,i,j),d).. abs(v(d,q) - v(d,i)*v(d,j)) =l= (2-wl(q,i)-wr(q,j))*B;
* Denominator
Termlimb.. sum(m, zb(m)) =l= BigM;
Coefflimb(m).. abs(be(m)) =l= zb(m)*B;
decompLb(q).. sum(Qi(q,i),gl(q,i)) =e= 1;
decompRb(q).. sum(Qj(q,j),gr(q,j)) =e= 1; 
valuesb(Qm(q,i,j),d).. abs(vb(d,q) - vb(d,i)*vb(d,j)) =l= (2-gl(q,i)-gr(q,j))*B;
* Final rational function surrogate
surrogate(d).. (a0+sum(m,a(m)*v(d,m)))/(b0+sum(m,be(m)*vb(d,m))) =e= ans(d);   

Model equ /all/;
Option minlp=BARON;
Solve equ using minlp minimizing er;
option decimals =8; display x,v.l,a.l,be.l,b0.l,a0.l,wl.l,wr.l,gl.l,gr.l,Y;



