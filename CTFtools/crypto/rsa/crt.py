from sympy.ntheory.modular import crt
import rsaATTACK
import gmpy2
def crtAtk(n_arr,c_arr,e):
    M = crt(n_arr,c_arr)
    rsaATTACK.send_res(gmpy2.iroot(M[0],e)[0])