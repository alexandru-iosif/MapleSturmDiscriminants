#####
with(PolynomialIdeals):
with(Groebner):
with(Student[MultivariateCalculus]):
with(Student[LinearAlgebra]):
with(combinat):

#####
SturmDiscriminants := module()
description "Sturm Discriminants";
#Author: Alexandru Iosif
option package;


#####
export SturmSequence, SturmDiscriminant, MonomialExponent, areAlgebraicallyIndependent, GenericPolynomial;


#####Main Functions:

SturmSequence;
SturmSequence := proc(f,x)
local c, n, i, fd, ci;
n := degree(f,x);
c := [];
if n >= 0 then
   c := [f];
   if n >= 1 then
      fd := diff(f,x);
      c := [f,fd];
      for i from 2 to n do
      	  ci := -rem(c[i-1],c[i],x);
      	  c := [op(c),ci];
      end do;
   end if;
end if;
return c;
end proc;

SturmDiscriminant;
SturmDiscriminant := proc(J,variablesf)
local c, t, Jt, parametersf, gensJt, SD, SDt, StSeq, coeffJt, p, d, i;
parametersf := IdealInfo[Variables](J) minus variablesf;
c := [];
SD := [];
for t in variablesf do
    Jt := EliminationIdeal(J,{op(parametersf),t});
    gensJt := IdealInfo[Generators](Jt);
    if numelems(gensJt) > 1 then
       return "Error: One of the elimination ideals is not principal. Please use another method to compute the discriminant.";
    end if;
    coeffJt := {coeffs(op(gensJt),t)};
    d := degree(op(gensJt),t);
    if areAlgebraicallyIndependent(coeffJt,parametersf) = true and d + 1 = numelems(coeffJt) then
       p := GenericPolynomial(d,t);
       StSeq := {op(SturmSequence(p,t))};
       SDt := {seq(coeffs(collect(numer(StSeq[i]),t,'distributed'),t),i=1..numelems(numer(StSeq)))}  union {seq(coeffs(denom(StSeq)[i],t),i=1..numelems(denom(StSeq)))} minus {1};
       SDt := [[op(SDt)],[seq({coeffs(p,t)}[i] = coeffJt[i],i=1..d+1)]]
#      The following is very slow. I comment it out;       
#      StSeq := subs({seq({coeffs(p,t)}[i] = coeffJt[i],i=1..d+1)},StSeq);
    else
       StSeq := {op(SD),op(SturmSequence(op(gensJt),t))};
       SDt := {seq(coeffs(collect(numer(StSeq)[i],t,'distributed'),t),i=1..numelems(numer(StSeq)))}  union {seq(coeffs(denom(StSeq)[i],t),i=1..numelems(denom(StSeq)))};
       SDt := [[op(SDt)],[]];
    end if;
#    SDt := {seq(coeffs(numer(StSeq)[i],t),i=1..numelems(numer(StSeq)))}  union {seq(coeffs(denom(StSeq)[i],t),i=1..numelems(denom(StSeq)))};
    SD := [op(SD),SDt];
end do;
return SD;
end proc;

MonomialExponent;
MonomialExponent := proc(m,variablesm)
local variablesmnew, c, variablesProduct, n, mu, i, exponentmu;
variablesmnew := [op(variablesm)];
n := numelems(variablesmnew);
if n = 1 then
   return [degree(m)];
end if;
exponentmu :=[];
variablesProduct := product(variablesmnew[i], i=1..n);
mu := sort(m*variablesProduct^2,variablesmnew,plex);
c := [op(mu)];
for i from 1 to n do
    exponentmu := [op(exponentmu), [op(c[i])][2]];
end do;
return exponentmu - [seq(2, i = 1 .. n)];
end proc;

areAlgebraicallyIndependent;
areAlgebraicallyIndependent := proc(L,variablesL)
local n,l;
n := numelems(L);
if n > numelems(variablesL) then
   return false;
end if;
for l in chose([op(variablesL)],n) do
    if Jacobian([op(L)],l) <> 0 then
       return true;
    end if;
end do;
return false;
end proc;

GenericPolynomial;
GenericPolynomial := proc(n,x)
local p, a, i;
a = {seq(a[i],i=0..n)};
p := sum(a[i]*x^i,i=0..n);
return p;
end proc;

end module; # SturmDiscriminants


#####Bibliography:
#[I19] A. Iosif, Algebraic methods for the study of multistationarity in mass-action networks, PhD thesis, Otto-von-Guericke-Universit\"at Magdeburg (2019)