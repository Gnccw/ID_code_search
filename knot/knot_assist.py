from cao_pylib import cipher_sbox
from cao_pylib import sbox
from cao_pylib import sbox_model
from cao_pylib import operation
import sys
sys.path.append('c:\\Users\\Gnccw\\Desktop\\imposs difference')
import main

'''
knot_sbox=cipher_sbox.knot_sbox
posspoint=main.tri_pset2bin_pset(main.undisturb_sbox(knot_sbox))
print(posspoint)
#main.test_ine('knot\\sbox_ine_from_sage.txt',posspoint)
'''
fp=open('knot/all_imposs_diff.txt','r')
s=fp.read()
fp.close()

inputdiffs=[]
s=s.split('\n')[:len(s.split('\n'))-1]
for v in s:
    if v[0]=='i':
        inputdiffs.append(v[10:])

for inputdiff in inputdiffs:
    res=''
    for j in range(64):
        tmp=[]
        for i in range(4):
            tmp.append(eval(inputdiff[(3-i)*64+j]))
        tmpnum=operation.binlist2num(tmp)
        res+=str(tmpnum)
    print(res)

            

        
