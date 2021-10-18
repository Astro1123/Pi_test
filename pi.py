#import numpy as np
import mpmath
import time

def fact(x):
	y=1
	for i in range(x):
		y*=(i+1)
	return y

def Newton(x):
    return x-mpmath.sin(x)/mpmath.cos(x)

def Newton_subroutine(x,cntmax):
	if cntmax == -1:
		while True:
			y=Newton(x)
			if x == y:
				x = y
				break
			x = y
	else:
		for i in range(cntmax):
			y=Newton(x)
			if x == y:
				x = y
				break
			x = y
	return x

def regular_polygon(ll,n):
      sl=ll/mpmath.sqrt(1.0+pow(ll/pow(2.0,n+1),2))
      return 2.0/(1.0/sl+1.0/ll)

def regular_polygon_subroutine(x,cntmax):
	i=1
	if cntmax == -1:
		while True:
			y=regular_polygon(x,i)
			if x == y:
				x = y
				i+=1
				break
			x = y
			i+=1
	else:
		for j in range(cntmax):
			y=regular_polygon(x,i)
			if x == y:
				x = y
				i+=1
				break
			x = y
			i+=1
	i-=1
	y=x/mpmath.sqrt(1.0+pow(x/pow(2.0,i+2),2))
	return x, y

def Accelerated(cntmax):
	x=0
	if cntmax == -1:
		i=0
		y=0
		while True:
			fact=1
			for j in range(i, 0, -1):
				fact=fact*(j*j/(2+2*j)/(1+2*j))
			x+=fact
			z = y
			y = mpmath.sqrt(x*9)
			i+=1
			if y == z:
				break
	else:
		for i in range(cntmax):
			fact=1
			for j in range(i, 0, -1):
				fact=fact*(j*j/(2+2*j)/(1+2*j))
			x+=fact
		y = mpmath.sqrt(x*9)
	return y

	
# ----------------------------- Main -----------------------------

print('precision : ', end='')
try:
	mpmath.mp.dps = int(input())
except ValueError:
	mpmath.mp.dps = 100

tmp = mpmath.mp.dps

print('')
print('count : ', end='')
try:
	cntmax = int(input())
except ValueError:
	cntmax = 10

print('')
print('input : ')
print('\tprecision\t : ',mpmath.mp.dps)
print('\tcount\t\t : ',cntmax)

print('')
print('Pi : ')
pi = mpmath.acos(-1)
print(pi)

print('')
print('Newton method : ')

stime = time.time()
x=Newton_subroutine(3,cntmax)
etime = time.time()

print(x)
print('( time : ',etime-stime,' sec )')
mpmath.mp.dps = 10
print('(error : ',x-pi,' )')
mpmath.mp.dps = tmp
print('')
print('Polygon : ')

stime = time.time()
ll,sl=regular_polygon_subroutine(4,cntmax)
etime = time.time()

print('sl : ')
print(sl)
print('ll : ')
print(ll)
print('( time : ',etime-stime,' sec )')
mpmath.mp.dps = 10
print('(error ( sl ) : ',sl-pi,' )')
print('(error ( ll ) : ',ll-pi,' )')
mpmath.mp.dps = tmp
print('')
print('Acceleration method : ')

stime = time.time()
x=Accelerated(cntmax)
etime = time.time()
print(x)
print('( time : ',etime-stime,' sec )')
mpmath.mp.dps = 10
print('(error : ',x-pi,' )')
mpmath.mp.dps = tmp
print('')