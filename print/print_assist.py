from cao_pylib import sbox_model
from cao_pylib import cipher_sbox
cnf="(0'+6)(2'+8)(4'+10)(5'+6)(1'+8)(3'+10)(3'+6)(5'+8)(1'+10)(0'+1)(2'+3)(1+3+4+10')(4'+5)(1+2+5+8')(0+3+5+6')(10'+11)(7+9')(9+11')(8'+11)(7'+8+10)"
cnf_inv="(0'+6)(2'+8)(4'+10)(5'+6)(1'+8)(3'+10)(3'+6)(5'+8)(1'+10)(0'+1)(2'+3)(1+3+4+10')(4'+5)(1+2+5+8')(0+3+5+6')(10'+11)(7+9')(9+11')(8'+11)(7'+8+10)"

print_sbox=cipher_sbox.print_sbox
p=[[0,1,3,6,7,4,5,2],
   [0,1,7,4,3,6,5,2],
   [0,3,1,6,7,5,4,2],
   [0,7,3,5,1,4,6,2]]

s0=[0, 1, 6, 5, 2, 7, 4, 3]
s1=[0, 1, 2, 7, 6, 5, 4, 3]
s2=[0, 6, 1, 5, 2, 4, 7, 3]
s3=[0, 2, 6, 4, 1, 7, 5, 3]

sbox_model.gen_ine_from_cnf(cnf,'print/print_sbox_ine.txt')
sbox_model.gen_ine_from_cnf(cnf_inv,'print/print_sbox_inv_ine.txt')

def pp(x):
    res=[0 for i in range(8)]
    for i in range(8):
        res[i]=print_sbox[x[i]]
    return res

print(pp(p[3]))