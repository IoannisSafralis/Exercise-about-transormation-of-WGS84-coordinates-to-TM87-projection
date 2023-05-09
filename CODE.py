import math as mt
import numpy as np
import statistics as st

def deg2rad(angle):
	return angle*mt.pi/180

def rad2deg(angle):
	return angle*180/mt.pi

X19 =4474197.366 
Y19=1897206.453 
Z19=4117253.641 
X30 =4464516.326 
Y30=1895649.704 
Z30=4128390.466 
X33 =4456830.186 
Y33=1912485.621 
Z33=4128945.737 
X35 =4467446.206 
Y35=1878147.197 
Z35=4132877.735 
X59 =4459394.586 
Y59=1887468.477 
Z59=4137485.151

#SYNISTOSES METATHESHS APO WGS84 -> EGSA87
dx = 199.723*np.ones(5)
dy = -74.030*np.ones(5)
dz = -246.018*np.ones(5)

X = [X19,X30,X33,X35,X59] + dx
Y = [Y19,Y30,Y33,Y35,Y59] + dy
Z = [Z19,Z30,Z33,Z35,Z59] + dz
i_d = [19,30,33,35,59] 

#WGS84
a= 6378137
e2 = 0.006694380
e_2 = 0.006739497

legsa=[]
fegsa=[]
hegsa=[]
Xegsa=[]
Yegsa=[]
Zegsa=[]
Np = []

for i in range(0,len(i_d)):
	print('Σημείο',i_d[i])
	P = mt.sqrt(X[i]**2+Y[i]**2)

	ln = rad2deg(mt.atan(Y[i]/X[i]))
	#****************************************************************
	fp = rad2deg(mt.atan((Z[i]*(1+e_2)/P)))
	#print(fp)

	w = mt.sqrt(1-e2*(mt.sin(deg2rad(fp)))**2)
	N = a/w
	fp1 = rad2deg(mt.atan((Z[i]+e2*N*mt.sin(deg2rad(fp)))/P))
	#print(fp1)

	w = mt.sqrt(1-e2*(mt.sin(deg2rad(fp1)))**2)
	N = a/w
	fp2 = rad2deg(mt.atan((Z[i]+e2*N*mt.sin(deg2rad(fp1)))/P))
	
	#print(fp2)
	#****************************************************************

	w = mt.sqrt(1-e2*(mt.sin(deg2rad(fp2)))**2)
	N = a/w
	fn = rad2deg(mt.atan((Z[i]+e2*N*mt.sin(deg2rad(fp2)))/P))
	print('φ =',fn)
	print('λ =',ln)

	h = Z[i]/mt.sin(deg2rad(fn)) - (1-e2)*N
	print('h =',h)
	print('*******************')

	fegsa.append(fn)
	legsa.append(ln)
	hegsa.append(h)

	Np.append(N)

	Xeg = (N+h)*mt.cos(deg2rad(fn))*mt.cos(deg2rad(ln))
	Yeg = (N+h)*mt.cos(deg2rad(fn))*mt.sin(deg2rad(ln))
	Zeg = ((1-e2)*N+h)*mt.sin(deg2rad(fn))

	Xegsa.append(Xeg)
	Yegsa.append(Yeg)
	Zegsa.append(Zeg)

#TM87
m0 = 0.9996
lo = deg2rad(24)
fo = 0
No = 0 
Eo = 500000
Epros = []
Npros = []

for i in range(0,len(i_d)):
	print('Σημείο',i_d[i])
	t = mt.tan(deg2rad(fegsa[i]))
	h2 = e_2*(mt.cos(deg2rad(fegsa[i])))**2
	dl = deg2rad(legsa[i]) - lo

	A0 = 1-1/4*e2-3/64*e2**2-5/256*e2**3-175/16384*e2**4
	A2 = 3/8*e2*(1+1/4*e2+15/128*e2**2+35/512*e2**3)
	A4 = 15/256*e2**2*(1+3/4*e2+35/64*e2**2)
	A6 = 35/3072*e2**3*(1+5/4*e2)
	Sf = a*(A0*deg2rad(fegsa[i])-A2*mt.sin(2*deg2rad(fegsa[i]))+A4*mt.sin(4*deg2rad(fegsa[i]))-A6*mt.sin(6*deg2rad(fegsa[i])))

	T1 = m0*Sf
	T2 = mt.sin(deg2rad(fegsa[i]))*mt.cos(deg2rad(fegsa[i]))/2
	T3 = mt.sin(deg2rad(fegsa[i]))*mt.cos(deg2rad(fegsa[i]))**3/24*(5-t**2+9*h2+4*h2**2)
	T4 = mt.sin(deg2rad(fegsa[i]))*(mt.cos(deg2rad(fegsa[i])))**5*(61-58*t**2+t**4+270*h2-330*t**2*h2+445*h2**2+324*h2**3-680*t**2*h2**2+88*h2**4-660*t**2*h2**3-192*t**2*h2**4)/720
	T5 = mt.sin(deg2rad(fegsa[i]))*(mt.cos(deg2rad(fegsa[i])))**7*(1385-3111*t**2+543*t**4-t**6)/40320
	Nn = No + T1 + m0*Np[i]*(T2*(dl)**2 + T3*(dl)**4 + T4*(dl)**6 + T5*(dl)**8)

	T6 = mt.cos(deg2rad(fegsa[i]))
	T7 = (mt.cos(deg2rad(fegsa[i])))**3*(1-t**2+h2)/6
	T8 = (mt.cos(deg2rad(fegsa[i])))**5*(5-18*t**2+t**4+14*h2-58*t**2*h2+13*h2**2+4*h2**3-64*t**2*h2**2-24*t**2*h2**3)/120
	T9 = (mt.cos(deg2rad(fegsa[i])))**7*(61-479*t**2+179*t**4-t**6)/5040
	En = Eo + m0*Np[i]*(T6*dl + T7*(dl)**3 + T8*(dl)**5 + T9*(dl)**7)
	
	Epros.append(En)
	Npros.append(Nn)

	print('E° =',En)
	print('N° =',Nn)
	print('*******************')

Est = [413246.537,415761.461,434262.815,398581.144,410371.156]
Nst = [4479171.649,4493790.178,4494324.815,4500086.407,4505901.944]

mEn = st.mean(Epros)

mNn = st.mean(Npros)

mEnn=[]
mEnn = mEn*np.ones((len(Epros)))

mNnn=[]
mNnn = mNn*np.ones((len(Npros)))

xi_ = Epros - mEnn
yi_ = Npros	- mNnn

P1 = sum(np.multiply(xi_, Est))
P2 = sum(np.multiply(yi_, Nst))
P3 = sum(np.multiply(yi_, Est))
P4 = sum(np.multiply(xi_, Nst))
P = sum(np.square(xi_)+np.square(yi_))

c = (P1 + P2)/P
d = (P3 - P4)/P

sx = sum(Est)/len(Est)
tx = sx - c*mEn - d*mNn
sy = sum(Nst)/len(Nst)
ty = sy + d*mEn - c*mNn

print('Παράμετροι μετασχηματισμού ομοιότητας:')
print('c =',c)
print('d =',d)
print('tx =',round(tx,3),'m')
print('ty =',round(ty,3),'m')

s = np.array([sx*np.ones(len(Est)),sy*np.ones(len(Nst))])
x_y_ = np.array([xi_,yi_])
cd = [[c,d],[-d,c]]
xy = [Est,Nst]

u = np.dot(cd,x_y_)
v = xy - s - u 

vx, vy = np.vsplit(v, 2)
#print(vx)
#print(vy)

rms = np.sqrt(np.mean(np.square(vx)+np.square(vy)))
print('RMS =',round(rms,4))

ETR = Est - vx
NTR = Nst - vy
ENTR = np.array([ETR,NTR])

Emt = [425668.532,423598.756,409679.172,410282.656,419271.864,427174.427,401065.602,413679.389]
Nmt = [4472158.023,4480091.718,4473286.345,4480023.339,4487912.070,4502251.541,4496687.333,4499605.765]


s = np.array([sx*np.ones(len(Emt)),sy*np.ones(len(Nmt))])
mEnn=[]
mEnn = mEn*np.ones((len(Emt)))

mNnn=[]
mNnn = mNn*np.ones((len(Nmt)))

ximt_ = Emt - mEnn
yimt_ = Nmt	- mNnn
xmt_ymt_ = np.array([ximt_,yimt_])

ENmTR = np.dot(cd,xmt_ymt_) + s

print(ENmTR)