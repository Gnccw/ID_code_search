from cao_pylib import undisturb
from cao_pylib import cipher_sbox
from ascon.ascon import *
from cao_pylib import main

fp=open('ascon_result.txt')
s=fp.read()
fp.close()

s=s.replace(' ', '')
s=s.replace('2', '?')

s=s.split('\n')

ss=''
for i in range(len(s)):
    s[i]="'"+s[i]+"'"+'\n'
    ss+=s[i]

fp=open('test.txt')
s=fp.read()
fp.close()

s=s.replace("'", '')
s=s.replace('?', '\s')
print(s)