Sets
*Total number of terms
m /1*5/
*Data point d
d 
*Total numer of input variables
N(m) 
*Total number of quadratic terms
q(m)
alias(m,i)
alias(m,j);
Sets
* Mapping valid ranges of i and j onto the corrosponding value of q
Qi(m,i) / 3. 1*2
          4. 1*3
          5. 1*4/
Qj(m,j) /set.Qi/
* Double mapping to produce every possible valid combination of i and j for the corrosponding value of q
Qm(m,i,j) /3.1*2.1*2
           4.1*3.1*3
           5.1*4.1*4/
* Valid combinations of q,i,j for various symmetry cuts           
QIi(m,i,j)  /3.1.2
             4.1.2*3
             4.2.3
             5.1.2*4
             5.2.3*4
             5.3.4/
Qi2(m,i) /3.1
          4.1*2
          5.1*3/;

Parameters
x(d,N) *Input data
Y(d)  *Response data
B * Coefficient limit
BigM *Maximum number of terms;

* Imports predictor data from csv;
$call csv2gdx predictor_data.csv id=x index=1 values=2..lastCol useheader=y
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
er Error
a0 Inital coefficient 
v(d,m) Value of terms;
Binary Variables
z(m) *Binary variable controlling each term
wl(m,i) *Left component binary variable 
wr(m,j) * Right component binary variable; 


B = 10000;
BigM = 5;
*Sets values of input predictor terms to the value of the corrosponding predictor at data point d
v.fx(d,N) = x(d,N);
*Initialises the program as a standard linear regression
z.l(q) = 0;
* Function to find
Y(d) = 7+2*x(d,'1')+4.1*x(d,'2')+x(d,'1')*x(d,'2')+sqr(x(d,'1'));
 
Equations
error,
Termlim,
Coefflim,
decompL,decompR,
values
Termsymmetry
RepeatTerm
Ordersymmetry;
* Least squares objective function
error.. sum(d, sqr(Y(d) - a0 - sum(m, a(m)*v(d,m)))) =e= er;
* Limits the total number of terms to BigM
Termlim.. sum(m, z(m)) =l= BigM;
* Forces coefficients to 0 for unselected terms and bounds them between +/- b otherwise
Coefflim(m).. abs(a(m)) =l= z(m)*B;
* Ensures a left and right hand operand are selected for monomial decomposition
decompL(q).. sum(Qi(q,i),wl(q,i)) =e= 1;
decompR(q).. sum(Qj(q,j),wr(q,j)) =e= 1;
* Enforces the identity of constructed terms as the product of terms selected by binary variables wl ,wr
values(Qm(q,i,j),d).. abs(v(d,q) - v(d,i)*v(d,j)) =l= (2-wl(q,i)-wr(q,j))*B;
* Helps prevent repeat terms
RepeatTerm(i,j).. sum(Qm(q,i,j),wl(q,i)*wr(q,j))=l=1;
* Helps remove intra-term sysmmetry
Termsymmetry(Qi2(q,i)).. wl(q,i) =l= 1 - (sum(QIi(q,i,j),wr(q,j)));
* Help remove inter-term symmetry by imposing an order on the terms
Ordersymmetry(q+1).. abs(v('1',q+1))-abs(v('1',q)) =g= 0 ;
    

Model equ /all/;
Option minlp=BARON;
Solve equ using minlp minimizing er;
display x,v.l,a.l,a0.l,wl.l,wr.l,Y;



