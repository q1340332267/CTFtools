import gmpy2
import os
import rsa
import requests
from math import log2
from Crypto.Util.number import *
from re import findall
from primefac import williams_pp1,pollard_pm1
def get_list_from_file(file,name,split=' = '):
    '''
    文件名
    所需的数组
    分隔符,默认 = 
    '''
    arr = []
    with open("{}".format(file),"rb")as f:
        res = findall("{} = [0-9]+".format(name),f.read().decode())
        for i in range(len(res)):
            arr.append(int(res[i].split(split)[1]))
        return arr
def send_res(*args):
    '''
    1. 传m 并打印
    2. 传c,d,n 打印
    '''
    if (len(args) == 1):
        print(long_to_bytes(args[0]).decode())
        exit
    if (len(args) == 3):
        print(long_to_bytes(gmpy2.powmod(args[0],args[1],args[2])).decode())
        exit

def factor_n(n):
    '''
    return array_p
    '''
    p = []
    res = requests.get("http://factordb.com/api?query={}".format(n))
    tmp = list(eval(res.text).items()) #eval将string转dict
    status = tmp[1][1]
    
    if (status == 'FF'):
        f_p = tmp[2:][0][1]
        for i in f_p:
            p.append(int(i[0]))
        return(p)
    else:
        exit("Can't factor,change method")
def p_plus_1(n):
    '''
    return array_p
    '''
    p = []
    tem = williams_pp1(n)
    p.append(tem)
    p.append(n//tem)
    return p
def p_minus_1(n):
    '''
    return array_p
    '''
    p = []
    tem = pollard_pm1(n)
    p.append(tem)
    p.append(n//tem)
    return p

def e_phi_ngcd(array_p,e,c):
    phi_N = 1
    n = 1
    for i in array_p:
        if (gmpy2.gcd(i - 1,e) == 1):
            phi_N *= (i-1)
            n *= i
    d = gmpy2.invert(e,phi_N)
    send_res(c,d,n)

def low_e(c,e,n):
    '''
    c,e,n
    '''
    k = 0
    while True:
        m,flag = gmpy2.iroot(c+k*n,e)
        k += 1
        if flag:
            send_res(m)
            break

class wiener:
    def transform(self,x,y):       #使用辗转相除将分数 x/y 转为连分数的形式
        res=[]
        while y:
            res.append(x//y)
            x,y=y,x%y
        return res
        
    def continued_fraction(self,sub_res):
        numerator,denominator=1,0
        for i in sub_res[::-1]:      #从sublist的后面往前循环
            denominator,numerator=numerator,i*numerator+denominator
        return denominator,numerator   #得到渐进分数的分母和分子，并返回

        
    #求解每个渐进分数
    def sub_fraction(self,x,y):
        res=self.transform(x,y)
        res=list(map(self.continued_fraction,(res[0:i] for i in range(1,len(res)))))  #将连分数的结果逐一截取以求渐进分数
        return res

    #以上是获得e/n的连分数

    def get_pq(self,a,b,c):      #由p+q和pq的值通过维达定理来求解p和q
        par=gmpy2.isqrt(b*b-4*a*c)   #由上述可得，开根号一定是整数，因为有解
        x1,x2=(-b+par)//(2*a),(-b-par)//(2*a)
        return x1,x2

    def wienerAttack(self,e,n):
        for (d,k) in self.sub_fraction(e,n):  #用一个for循环来注意试探e/n的连续函数的渐进分数，直到找到一个满足条件的渐进分数
            if k==0:                     #可能会出现连分数的第一个为0的情况，排除
                continue
            if (e*d-1)%k!=0:             #ed=1 (\pmod φ(n)) 因此如果找到了d的话，(ed-1)会整除φ(n),也就是存在k使得(e*d-1)//k=φ(n)
                continue
            
            phi=(e*d-1)//k               #这个结果就是 φ(n)
            px,qy=self.get_pq(1,n-phi+1,n)
            if px*qy==n:
                p,q=abs(int(px)),abs(int(qy))     #可能会得到两个负数，负负得正未尝不会出现
                d=gmpy2.invert(e,(p-1)*(q-1))     #求ed=1 (\pmod  φ(n))的结果，也就是e关于 φ(n)的乘法逆元d
                return d, p, q
        print("winer失效")
    pass
def winer_attack(*args):
    '''
    c,e,n
    d,n  print(p,q)
    '''
    if len(args) == 3:
        d = wiener().wienerAttack(args[1],args[2])
        send_res(args[0],d,args[2])
    if len(args) == 2:
        _,p,q = wiener().wienerAttack(args[0],args[1])
        print (p,q)

def n_ngcd(c_arr,e,n_arr):
    flag = False
    for i in range(len(n_arr)):
        for j in range(len(n_arr)):
            if (i != j):
                p = gmpy2.gcd(n_arr[i],n_arr[j])
                if p != 1:
                    q = n_arr[i] // p
                    d = gmpy2.invert(e,(p-1)*(q-1))
                    send_res(c_arr[i],d,n_arr[i])
                    flag = True
                    break
        if flag:
            break            


def dp_leak(c,e,n,dp):
    '''
    c,e,n,dp
    '''
    if e < 1000000 :
        for x in range(1,e):
            p = (dp*e - 1)//x +1
            if gmpy2.gcd(n,p) != 1:
                q = n // p
                d = gmpy2.invert(e,(p-1)*(q-1))
                send_res(c,d,n)
                break
    else:
        m = 1000000007
        p = gmpy2.gcd(gmpy2.gmpy2.powmod(m, e*dp, n) - m, n)
        if p != 1:
            q = n // p
            d = gmpy2.invert(e, (p - 1) * (q - 1))
            m = gmpy2.powmod(c, d, n)
            send_res(m)

def dpdq_leak(c,p,q,dp,dq):
    '''
    c,p,q,dp,dq
    '''
    invq = gmpy2.invert(p, q)
    m1 = gmpy2.powmod(c, dp, p)
    m2 = gmpy2.powmod(c, dq, q)
    m = (((m2 - m1)*invq) % q)*p + m1
    send_res(m)

def same_mod(e1,e2,n,c1,c2):
    _,s,t = gmpy2.gcdext(e1,e2)
    m = pow(c1,s,n)*pow(c2,t,n)%n
    send_res(m)

class rabin_Class:
    def rabin_decrypt(this, c, p, q, e):
        n = p*q
        x0=gmpy2.invert(p,q)
        x1=gmpy2.invert(q,p)
        cs = [c]
        for _ in range(int(log2(e))):
            ps = []
            for c2 in cs:
                r = pow(c2, (p + 1) // 4, p)
                s = pow(c2, (q + 1) // 4, q)

                x = (r * x1 * q + s * x0 * p) % n
                y = (r * x1 * q - s * x0 * p) % n
                if x not in ps:
                    ps.append(x)
                if n - x not in ps:
                    ps.append(n - x)
                if y not in ps:
                    ps.append(y)
                if n - y not in ps:
                    ps.append(n - y)
            cs = ps
        return cs
def rabin(c, p, q, e):
    rb = rabin_Class().rabin_decrypt(c, p, q, e)
    for i in range(4):
        try:
            send_res(rb[i])
        except:
            pass