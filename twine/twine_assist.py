from cao_pylib import sbox
from cao_pylib import sbox_model
F=[[0,0],[1,2],[2,2],[3,4],[4,4]]
xor=[[0,0,0],
     [0,1,1],
     [0,2,2],
     [0,3,3],
     [0,4,4],
     [1,0,1],
     [1,1,1],
     [1,1,0],
     [1,2,3],
     [1,3,4],
     [1,3,2],
     [1,4,4],
     [2,0,2],
     [2,1,3],
     [2,2,4],
     [2,3,4],
     [2,4,4],
     [3,0,3],
     [3,1,4],
     [3,1,2],
     [3,2,4],
     [3,3,4],
     [3,4,4],
     [4,0,4],
     [4,1,4],
     [4,2,4],
     [4,3,4],
     [4,4,4]]

bin_F=[]
bin_xor=[]

for i in F:
    tmp=[]
    for j in i:
        if j==0:
            tmp+=[0,0,0]
        if j==1:
            tmp+=[0,0,1]
        if j==2:
            tmp+=[0,1,0]
        if j==3:
            tmp+=[0,1,1]
        if j==4:
            tmp+=[1,0,0]
    bin_F.append(tmp)

for i in xor:
    tmp=[]
    for j in i:
        if j==0:
            tmp+=[0,0,0]
        if j==1:
            tmp+=[0,0,1]
        if j==2:
            tmp+=[0,1,0]
        if j==3:
            tmp+=[0,1,1]
        if j==4:
            tmp+=[1,0,0]
    bin_xor.append(tmp)

sbox.write_patterns_csv(bin_F,'twine//bin_F.csv')
sbox.write_patterns_csv(bin_xor,'twine//bin_xor.csv')

cnf_F="(5')(0'+3)(1+2+4')(3'+4')(1+2'+4)(1'+3+4)(0'+1')(0+2+3')(1'+2'+4')"
cnf_xor="(3'+4')(0'+1')(6'+7')(3'+5')(0'+2')(3'+6)(0'+6)(2+5+8')(1+4+7')(4'+6+7)(1'+6+7)(6'+8')(2+5'+6+8)(2'+5+6+8)(0+1+3+4+6')(1'+4'+7')(0+1+4'+5+6')(1'+3+4+5+6')(2'+5'+7'+8')(0+1+2+5'+6')(1'+2+4+5'+6')"
sbox_model.gen_ine_from_cnf(cnf_F,'twine//ine//F_ine.txt')
sbox_model.gen_ine_from_cnf(cnf_xor,'twine//ine//xor_ine.txt')