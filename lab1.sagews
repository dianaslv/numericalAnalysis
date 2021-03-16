
#Problema1
#MacLauren -> Taylor in jurul pct 0

var('x a') #declar varibalilele simbolice
assume(a > 0)
f(x) = sqrt(x+a)
show(f.taylor(x,0,10))


#Problema2

var('n')
R_n = 3/factorial(n+1)
for order in range(1,200):
    if R_n(order) < 10^(-3):
        print(str(order))
        break

#Problema3

var('n')
R_n(n)=binomial(1/3,n)*10^(-3*n)
for order in range(1,2000):
    if abs(R_n(order))<10^(-12):
        print(str(order))
        break

#Problema4

#METODA1: integram (seria taylor pt exp_t)
var('t')
assume(t>0)
exp_t(t) = exp(-t^2)
exp_taylor(t) = exp_t.taylor(t,0,5)
erf1(x) = 2/(sqrt(pi))*exp_taylor(t).integral(t,0,x)
print(erf1)
print(str(RR(erf1(1))))

#METODA2: calc seria talor pt erf in jurul lui 0
erf2(x) = erf(x).taylor(x,0,5)
print(erf)
print(erf2)
print(str(RR(erf2(1))))

#Problema5

#Deduceti seria Taylor pentru ln(1+x)
taylor(ln(1+x),x,0,21)  #afiseaza seria taylor (elementele)


#aproximati ln 2 folosind primii 8 termeni ai seriei Taylor
var('x')
assume(x>0)
ln_taylor(x)=ln(1+x).taylor(x,0,8) #calculeaza seria taylor pt x termeni dati
print(str(RR(ln_taylor(1))))


#Cati termeni sunt necesari pentru a obtine ln 2 cu 5 zecimale corecte?
var('n')
assume(n>0)
R_n(n) = 1/(n+1)
order = 1
while True:
    if(abs(R_n(order))<10^(-5)):
        print(str(order))
        break
    order +=1

#idem pt ln(1+x)/(1-x)
taylor(ln((1+x)/(1-x)),x,0,21)


var('x')
assume(x>0)
ln_taylor(x)=ln((1+x)/(1-x)).taylor(x,0,8) #calculeaza seria taylor pt x termeni dati
print(str(RR(ln_taylor(1/3))))

var('x')
assume(x>0)
R_n(n) = 1/(3^n)*(2*(n+1))
order = 1
while True:
    if(abs(R_n(order))<10^(-5)):
        print(str(order))
        break
    order +=1

#Problema6
taylor(arctan(x),x,0,21)


#Cati termeni sunt necesari pentru a obtine π/4 cu 5 zecimale corecte.
var('x')
assume(x>0)
R_n(n) = 1/(2*(n+1))
order = 1
while True:
    if(abs(R_n(order))<10^(-5)):
        print(str(order))
        break
    order +=1

#Problema 7
var('x')
def macLuarinTrunch(f,order):
    return sum([(derivative(f,x,i)(x=0)/factorial(i))*x^i for i in range(order+1)])

T2f(x) = macLuarinTrunch(exp(x),2)
T3f(x) = macLuarinTrunch(exp(x),3)
T4f(x) = macLuarinTrunch(exp(x),4)
T5f(x) = macLuarinTrunch(exp(x),5)

p1=plot(exp(x),(x,-2,2))
p2=plot(T2f,(x,-2,2))
p3=plot(T3f,(x,-2,2))
p4=plot(T4f,(x,-2,2))
p5=plot(T5f,(x,-2,2))
p1+p2+p3+p4+p5


T2f(x) = macLuarinTrunch(ln(1+x),2)
T3f(x) = macLuarinTrunch(ln(1+x),3)
T4f(x) = macLuarinTrunch(ln(1+x),4)
T5f(x) = macLuarinTrunch(ln(1+x),5)

p1=plot(ln(1+x),(x,-0.99,2))
p2=plot(T2f,(x,-0.99,2))
p3=plot(T3f,(x,-0.99,2))
p4=plot(T4f,(x,-0.99,2))
p5=plot(T5f,(x,-0.99,2))
p1+p2+p3+p4+p5

#Problema 8
#Pade
def pade(f,m,k):
    c=[derivative(f,x,i)(x=0)/factorial(i) for i in range(m+k+1)]
    C = matrix.toeplitz([c[m+1] for i in range(k)], [c[m-i] if m-i>0 else 0 for i in range(1,k)])
    R = vector([-c[i] for i in range(m+1,m+k+1)])
    B= C.solve_right(R)
    B_coeff = [1] + [elem for elem in B]
    if m>k:
        for i in range(k+1,m+1):
            B_coeff.append(0)
    B = vector(B_coeff)
    A = [sum([c[j-l]*B[l] for l in range(j+1)]) for j in range(m+1)]
    return sum([A[i]*x^i for i in range(m+1)])/sum(B[i]*x^i for i in range(k+1))


R11(x) = pade(exp(x),1,1)
R22(x) = pade(exp(x),2,2)

p1 = plot(exp(x),(x,-1,1));
p2 = plot(R11(x),(x,-1,1));
p3 = plot(R22(x),(x,-1,1));
p1 +p2+p3









