from gurobipy import *
from cao_pylib import sbox_model
from cao_pylib import operation

def sbox_xor(m,x,y,z):
    for i in range(len(x)):
        tmp=x[i]+y[i]+z[i]
        sbox_model.gen_model_from_ine(m,tmp,'ine//xor_ine.txt')

def print_table(x):
    for xx in x:
        print(xx)

def sbox_3xor(m,x,y,z,w):
    for i in range(len(x)):
        tmp=x[i]+y[i]+z[i]+w[i]
        sbox_model.gen_model_from_ine(m,tmp,'ine//3xor_ine.txt')

def sbox_equ(m,x,y):
    for i in range(len(x)):
        m.addConstr(x[i][0]==y[i][0])
        m.addConstr(x[i][1]==y[i][1])

def list2value(x):
    if x==[0,0]:
        return 0
    if x==[0,1]:
        return 1
    if x==[1,1]:
        return 2
    else:
        return 'err'

def _get_value(skinny):
    xnum=len(skinny.x)
    snum=len(skinny.sout)
    rounds=skinny.r
    flag=skinny.flag
    sbox_size=skinny.sbox_size
    x=[0 for i in range(xnum)]
    sout=[0 for i in range(snum)]

    for v in skinny.m.getVars():
        if v.varName[0]=='x':
            x[eval(v.varName[2:len(v.varName)-1])]=int(v.x)
        if v.varName[0]=='s':
            sout[eval(v.varName[5:len(v.varName)-1])]=int(v.x)
    
    for r in range(rounds):
        tmpx=[[0 for i in range(4)] for j in range(4)]
        tmps=[[0 for i in range(4)] for j in range(4)]
        for row in range(4):
            for col in range(4):
                s=''
                for i in range(sbox_size):
                    s+=str(list2value([x[(r*128+(row*4+col)*sbox_size+i)*2],x[(r*128+(row*4+col)*sbox_size+i)*2+1]]))
                tmpx[row][col]=s
                s=''
                for i in range(sbox_size):
                    s+=str(list2value([sout[(r*128+(row*4+col)*sbox_size+i)*2],sout[(r*128+(row*4+col)*sbox_size+i)*2+1]]))
                tmps[row][col]=s
        print_table(tmpx)
        print('Sbox')
        print_table(tmps)
        print('linear')
    if flag==0:
        tmpx=[[0 for i in range(4)] for j in range(4)]
        for row in range(4):
            for col in range(4):
                s=''
                for i in range(sbox_size):
                    s+=str(list2value([x[(rounds*128+(row*4+col)*sbox_size+i)*2],x[(rounds*128+(row*4+col)*sbox_size+i)*2+1]]))
                tmpx[row][col]=s
        print_table(tmpx)


                

class Skinny():
    def __init__(self,rounds,blocksize,flag,index,value):
        self.r=rounds
        self.blocksize=blocksize
        self.index=index
        self.flag=flag
        if value==0:
            self.value=[0,0]
        if value==1:
            self.value=[0,1]
        self.sbox_size=int(self.blocksize/16)
        
        self.m=Model()
        self.x=self.m.addVars([i for i in range(self.blocksize*2*(self.r+1))],vtype=GRB.BINARY,name='x')
        self.sout=self.m.addVars([i for i in range(self.blocksize*2*(self.r))],vtype=GRB.BINARY,name='sout')
        self.a=self.m.addVars([i for i in range(self.r*16)],vtype=GRB.BINARY,name='at')

        constr=LinExpr()
        for i in range(self.blocksize*2):
            constr+=self.x[i]
        self.m.addConstr(constr>=1)
    
    def gen_constr(self):
        for r in range(self.r):
            xindex_start=self.blocksize*r*2

            xstate=[[0 for i in range(4)] for i in range(4)]
            sout_state=[[0 for i in range(4)] for i in range(4)]
            x1state=[[0 for i in range(4)]  for i in range(4)]
            a=[[0 for i in range(4)]  for i in range(4)]

            for row in range(4):
                for col in range(4):
                    
                    s=row*4+col
                    tmpx=[]
                    tmps=[]
                    tmpx1=[]
                    for i in range(self.sbox_size):
                        tmpx.append([self.x[xindex_start+(s*self.sbox_size+i)*2],self.x[xindex_start+(s*self.sbox_size+i)*2+1]])
                        tmps.append([self.sout[xindex_start+(s*self.sbox_size+i)*2], self.sout[xindex_start+(s*self.sbox_size+i)*2+1]])
                        tmpx1.append([self.x[xindex_start+(s*self.sbox_size+i+self.blocksize)*2],self.x[xindex_start+(s*self.sbox_size+i+self.blocksize)*2+1]])
                    self.m.update()
                    xstate[row][col]=tmpx
                    sout_state[row][col]=tmps
                    x1state[row][col]=tmpx1
                    if self.sbox_size==4:
                        tmp=[]
                        for i in range(4):
                            tmp=tmp+[tmpx[i][0],tmpx[i][1]]
                        for i in range(4):
                            tmp=tmp+[tmps[i][0],tmps[i][1]]
                        sbox_model.gen_model_from_ine(self.m,tmp,'skinny//skinny_sbox_4.txt')
                    if self.sbox_size==8:
                        for i in range(8):
                            for j in range(2):
                                tmp=[]
                                for k in range(8):
                                    tmp=tmp+[tmpx[k][0],tmpx[k][1]]
                                tmp+=[tmps[i][j]]
                                pos=i*2+j
                                sbox_model.gen_model_from_ine(self.m,tmp,'skinny//skinny'+str(pos)+'.txt')
            #线性层
            s1=[]
            for i in range(4):
                tmp=operation.roate_right(sout_state[i],i)
                s1.append(tmp)
            for col in range(4):
                sbox_3xor(self.m,s1[0][col],s1[2][col],s1[3][col],x1state[0][col])
                sbox_equ(self.m,s1[0][col],x1state[1][col])
                sbox_xor(self.m,s1[1][col],s1[2][col],x1state[2][col])
                sbox_xor(self.m,s1[0][col],s1[2][col],x1state[3][col])
        
        if self.flag==0:
            index=(self.r*self.blocksize+self.index)*2
            self.m.addConstr(self.x[index]==self.value[0])
            self.m.addConstr(self.x[index+1]==self.value[1])
        if self.flag==1:
            index=((self.r-1)*self.blocksize+self.index)*2
            self.m.addConstr(self.sout[index]==self.value[0])
            self.m.addConstr(self.sout[index+1]==self.value[1])
    
        
                            

for i in range(122,128):
    print(i)
    test=Skinny(5,128,0,i,0)
    test.gen_constr()
    test.m.optimize()
    if test.m.status!=GRB.Status.INFEASIBLE: 
        _get_value(test)
        break
       