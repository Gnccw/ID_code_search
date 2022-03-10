from gurobipy import *
from cao_pylib import sbox_model
from cao_pylib import operation

def f0(m,x,y):
    x1=operation.roate_left(x,1)
    x2=operation.roate_left(x,2)
    x3=operation.roate_left(x,7)
    for i in range(len(x)):
        sbox_model.gen_model_from_ine(m,x1[i]+x2[i]+x3[i]+y[i],'ine//3xor_ine.txt')

def f1(m,x,y):
    x1=operation.roate_left(x,3)
    x2=operation.roate_left(x,4)
    x3=operation.roate_left(x,6)
    for i in range(len(x)):
        sbox_model.gen_model_from_ine(m,x1[i]+x2[i]+x3[i]+y[i],'ine//3xor_ine.txt')

def equ(m,x,y):
    for i in range(len(x)):
        m.addConstr(x[i][0]==y[i][0])
        m.addConstr(x[i][1]==y[i][1])

def mand0(m,a,c,y):
    ta=[]
    tb=[]
    ty=[]
    tc=[]
    for i in range(8):
        ta.append(a[7-i])
        tb.append([0,0])
        ty.append(y[7-i])
    for i in range(7):
        tc.append(c[6-i])
    for i in range(8):
        if i==0:
            sbox_model.gen_model_from_ine(m,ta[i]+tb[i]+ty[i],'ine//xor_ine.txt')
            sbox_model.gen_model_from_ine(m,ta[i]+tb[i]+tc[i],'ine//2and_ine.txt')
        elif i<7:
            sbox_model.gen_model_from_ine(m,ta[i]+tb[i]+tc[i-1]+ty[i],'ine//3xor_ine.txt')
            #print(a[i]+b[i]+c[i-1]+y[i])
            sbox_model.gen_model_from_ine(m,ta[i]+tb[i]+tc[i-1]+tc[i],'arx//进位_3_ine.txt')
        else:
            sbox_model.gen_model_from_ine(m,ta[i]+tb[i]+tc[i-1]+ty[i],'ine//3xor_ine.txt')


def mand(m,a,b,c,y):
    ta=[]
    tb=[]
    ty=[]
    tc=[]
    for i in range(8):
        ta.append(a[7-i])
        #tb.append([0,0])
        tb.append(b[7-i])
        ty.append(y[7-i])
    for i in range(7):
        tc.append(c[6-i])
    for i in range(8):
        if i==0:
            sbox_model.gen_model_from_ine(m,ta[i]+tb[i]+ty[i],'ine//xor_ine.txt')
            sbox_model.gen_model_from_ine(m,ta[i]+tb[i]+tc[i],'ine//2and_ine.txt')
        elif i<7:
            sbox_model.gen_model_from_ine(m,ta[i]+tb[i]+tc[i-1]+ty[i],'ine//3xor_ine.txt')
            sbox_model.gen_model_from_ine(m,ta[i]+tb[i]+tc[i-1]+tc[i],'arx//进位_3_ine.txt')
        else:
            sbox_model.gen_model_from_ine(m,ta[i]+tb[i]+tc[i-1]+ty[i],'ine//3xor_ine.txt')


class Hight():
    def __init__(self,rounds,block_size,index,value):
        self.r=rounds
        self.index=index
        self.block_size=block_size
        if value==0:
            self.value=[0,0]
        if value==1:
            self.value=[0,1]
        
        self.m=Model()
        self.x=self.m.addVars([i for i in range((self.r+1)*128)],vtype=GRB.BINARY,name='x')
        self.c=self.m.addVars([i for i in range((self.r)*56)],vtype=GRB.BINARY,name='c')
        self.xout=self.m.addVars([i for i in range((self.r)*64)],vtype=GRB.BINARY,name='xout')
        self.mout=self.m.addVars([i for i in range((self.r)*32)],vtype=GRB.BINARY,name='mout')

        #输入2的个数
        self.d=self.m.addVars([i for i in range(64)],vtype=GRB.BINARY,name='d')
        for i in range(64):
            self.m.addConstr(self.x[i*2]>=self.d[i])
            self.m.addConstr(self.x[i*2+1]>=self.d[i])
            self.m.addConstr(self.x[i*2]+self.x[i*2+1]-self.d[i]<=1)
        obj=LinExpr()
        for i in range(64):
            obj+=self.d[i]
        self.m.setObjective(obj,GRB.MAXIMIZE)

        constr=LinExpr()
        for i in range(128):
            constr+=self.x[i]
        self.m.addConstr(constr>=1)

    def gen_constr(self):
        for r in range(self.r):
            xindex=r*64
            x=[]
            c=[]
            x1=[]
            xout=[]
            mout=[]
            for i in range(8):
                tmpx=[]
                tmpx1=[]
                for j in range(8):
                    tmpx.append([self.x[(xindex+i*8+j)*2],self.x[(xindex+i*8+j)*2+1]])
                    tmpx1.append([self.x[(xindex+i*8+j+64)*2],self.x[(xindex+i*8+j+64)*2+1]])
                x.append(tmpx)
                x1.append(tmpx1)
            for i in range(4):
                tmpc=[]
                tmpxout=[]
                for j in range(7):
                    tmpc.append([self.c[(28*r+i*7+j)*2],self.c[(28*r+i*7+j)*2+1]])
                for j in range(8):
                    tmpxout.append([self.xout[(32*r+i*8+j)*2],self.xout[(32*r+i*8+j)*2+1]])
                c.append(tmpc)
                xout.append(tmpxout)
            
            for i in range(2):
                tmpmout=[]
                for j in range(8):
                    tmpmout.append([self.mout[(16*r+i*8+j)*2],self.mout[(16*r+i*8+j)*2+1]])
                mout.append(tmpmout)

            for i in range(4):
                equ(self.m,x1[i*2],x[i*2+1])
            f0(self.m,x[1],xout[0])
            f1(self.m,x[3],xout[1])
            f0(self.m,x[5],xout[2])
            f1(self.m,x[7],xout[3])

            mand0(self.m,xout[0],c[0],mout[0])
            mand(self.m,xout[1],x[2],c[1],x1[1])
            mand0(self.m,xout[2],c[2],mout[1])
            mand(self.m,xout[3],x[6],c[3],x1[5])

            for i in range(8):
                sbox_model.gen_model_from_ine(self.m,x[0][i]+mout[0][i]+x1[7][i],'ine//xor_ine.txt')
                sbox_model.gen_model_from_ine(self.m,x[4][i]+mout[1][i]+x1[3][i],'ine//xor_ine.txt')
        index=(64*self.r+self.index)*2
        self.m.addConstr(self.x[index]==self.value[0],name='output0')
        self.m.addConstr(self.x[index+1]==self.value[1],name='output1')


class Hight_inv():
    def __init__(self,rounds,block_size,index,value):
        self.r=rounds
        self.index=index
        self.block_size=block_size
        if value==0:
            self.value=[0,0]
        if value==1:
            self.value=[0,1]
        
        self.m=Model()
        self.x=self.m.addVars([i for i in range((self.r+1)*128)],vtype=GRB.BINARY,name='x')
        self.c=self.m.addVars([i for i in range((self.r)*56)],vtype=GRB.BINARY,name='c')
        self.xout=self.m.addVars([i for i in range((self.r)*64)],vtype=GRB.BINARY,name='xout')
        self.mout=self.m.addVars([i for i in range((self.r)*32)],vtype=GRB.BINARY,name='mout')   

        #输入2的个数
        self.d=self.m.addVars([i for i in range(64)],vtype=GRB.BINARY,name='d')
        for i in range(64):
            self.m.addConstr(self.x[i*2]>=self.d[i])
            self.m.addConstr(self.x[i*2+1]>=self.d[i])
            self.m.addConstr(self.x[i*2]+self.x[i*2+1]-self.d[i]<=1)
        obj=LinExpr()
        for i in range(64):
            obj+=self.d[i]
        self.m.setObjective(obj,GRB.MAXIMIZE)    

        constr=LinExpr()
        for i in range(128):
            constr+=self.x[i]
        self.m.addConstr(constr>=1)

    def gen_constr(self): 
        for r in range(self.r):
            xindex=r*64
            x0=[]
            x1=[]
            c=[]
            xout=[]
            mout=[]
            for i in range(8):
                tmpx=[]
                tmpx1=[]
                for j in range(8):
                    tmpx.append([self.x[(xindex+i*8+j)*2],self.x[(xindex+i*8+j)*2+1]])
                    tmpx1.append([self.x[(xindex+i*8+j+64)*2],self.x[(xindex+i*8+j+64)*2+1]])
                x0.append(tmpx)
                x1.append(tmpx1)
            
            for i in range(4):
                tmpc=[]
                tmpxout=[]
                for j in range(7):
                    tmpc.append([self.c[(28*r+i*7+j)*2],self.c[(28*r+i*7+j)*2+1]])
                for j in range(8):
                    tmpxout.append([self.xout[(32*r+i*8+j)*2],self.xout[(32*r+i*8+j)*2+1]])
                c.append(tmpc)
                xout.append(tmpxout)

            for i in range(2):
                tmpmout=[]
                for j in range(8):
                    tmpmout.append([self.mout[(16*r+i*8+j)*2],self.mout[(16*r+i*8+j)*2+1]])
                mout.append(tmpmout)
            
            for i in range(4):
                equ(self.m,x0[i*2],x1[i*2+1])
            
            f0(self.m,x1[1],xout[0])
            f1(self.m,x1[3],xout[1])
            f0(self.m,x1[5],xout[2])
            f1(self.m,x1[7],xout[3])

            mand0(self.m,xout[0],c[0],mout[0])
            mand(self.m,xout[1],x0[1],c[1],x1[2])
            mand0(self.m,xout[2],c[2],mout[1])
            mand(self.m,xout[3],x0[5],c[3],x1[6])

            for i in range(8):
                sbox_model.gen_model_from_ine(self.m,x0[7][i]+mout[0][i]+x1[0][i],'ine//xor_ine.txt')
                sbox_model.gen_model_from_ine(self.m,x0[3][i]+mout[1][i]+x1[4][i],'ine//xor_ine.txt')

        index=(64*self.r+self.index)*2
        self.m.addConstr(self.x[index]==self.value[0],name='output0')
        self.m.addConstr(self.x[index+1]==self.value[1],name='output1')           



def _get_value_b(cla):
    x=[0 for i in range(len(cla.x))]
    c=[0 for i in range(len(cla.c))]
    xout=[0 for i in range(len(cla.xout))]
    mout=[0 for i in range(len(cla.mout))]

    for v in cla.m.getVars():
        if v.varName[0]=='x' and v.varName[1]!='o':
            x[eval(v.varName[2:len(v.varName)-1])]=int(v.x)
        if v.varName[0]=='c':
            c[eval(v.varName[2:len(v.varName)-1])]=int(v.x)
        '''
        if v.varName[0]=='x' and v.varName[1]=='o':
            xout[eval(v.varName[2:len(v.varName)-4])]=int(v.x)
        if v.varName[0]=='m':
            mout[eval(v.varName[2:len(v.varName)-4])]=int(v.x)
        '''
    for r in range(cla.r,-1,-1):
        print(r)
        s=''
        for i in range(64):
            if i%8==0:
                s+=' '
            if [x[r*128+i*2],x[r*128+i*2+1]]==[0,0]:
                s+='0'
            if [x[r*128+i*2],x[r*128+i*2+1]]==[0,1]:
                s+='1'
            if [x[r*128+i*2],x[r*128+i*2+1]]==[1,1]:
                s+='2'
        s+='\n'
        '''
        if r<lea.r:
            for i in range(93):
                if i%31==0:
                    s+=' '
                if [c[r*186+i*2],c[r*186+i*2+1]]==[0,0]:
                    s+='0'
                if [c[r*186+i*2],c[r*186+i*2+1]]==[0,1]:
                    s+='1'
                if [c[r*186+i*2],c[r*186+i*2+1]]==[1,1]:
                    s+='2'
        '''
        print(s)
        print()       
       
def _get_value_f(cla):
    x=[0 for i in range(len(cla.x))]
    c=[0 for i in range(len(cla.c))]
    xout=[0 for i in range(len(cla.xout))]
    mout=[0 for i in range(len(cla.mout))]

    for v in cla.m.getVars():
        if v.varName[0]=='x' and v.varName[1]!='o':
            x[eval(v.varName[2:len(v.varName)-1])]=int(v.x)
        if v.varName[0]=='c':
            c[eval(v.varName[2:len(v.varName)-1])]=int(v.x)
        '''
        if v.varName[0]=='x' and v.varName[1]=='o':
            xout[eval(v.varName[2:len(v.varName)-4])]=int(v.x)
        if v.varName[0]=='m':
            mout[eval(v.varName[2:len(v.varName)-4])]=int(v.x)
        '''
    for r in range(cla.r+1):
        print(r)
        s=''
        for i in range(64):
            if i%8==0:
                s+=' '
            if [x[r*128+i*2],x[r*128+i*2+1]]==[0,0]:
                s+='0'
            if [x[r*128+i*2],x[r*128+i*2+1]]==[0,1]:
                s+='1'
            if [x[r*128+i*2],x[r*128+i*2+1]]==[1,1]:
                s+='2'
        s+='\n'
        
        print(s)
        print()       

def Hight_forward(r,x):
    input_x=[]
    test=Hight(r,64,1,1)
    test.gen_constr()
    test.m.update()
    test.m.remove(test.m.getConstrByName('output0'))
    test.m.remove(test.m.getConstrByName('output1'))
    for i in range(len(x)):
        if x[i]!=' ':
            if x[i]=='0':
                input_x+=[0,0]
            if x[i]=='1':
                input_x+=[0,1]
            if x[i]=='2':
                input_x+=[1,1]
    for i in range(128):
        test.m.addConstr(test.x[i]==input_x[i])
    test.m.Params.OutputFlag=0
    test.m.optimize()
    if test.m.status!=GRB.Status.INFEASIBLE:
        _get_value_f(test)

def Hight_backward(r,x):
    input_x=[]
    test=Hight_inv(r,64,1,1)
    test.gen_constr()
    test.m.update()
    test.m.remove(test.m.getConstrByName('output0'))
    test.m.remove(test.m.getConstrByName('output1'))
    for i in range(len(x)):
        if x[i]!=' ':
            if x[i]=='0':
                input_x+=[0,0]
            if x[i]=='1':
                input_x+=[0,1]
            if x[i]=='2':
                input_x+=[1,1]
    for i in range(128):
        test.m.addConstr(test.x[i]==input_x[i])
    test.m.Params.OutputFlag=0
    test.m.optimize()
    if test.m.status!=GRB.Status.INFEASIBLE:
        _get_value_b(test)

'''
for i in range(8,64):
    print(i)           
    test=Hight(8,64,i,1)
    test_inv=Hight_inv(8,64,i,0)

    test.gen_constr()
    test.m.Params.OutputFlag=0
    test.m.optimize()

    
    test_inv.gen_constr()
    test_inv.m.Params.OutputFlag=0
    test_inv.m.optimize()

    if test.m.status!=GRB.Status.INFEASIBLE and test_inv.m.status!=GRB.Status.INFEASIBLE:
        _get_value_f(test)
        _get_value_b(test_inv)
        print(i,'yes') 
        break
'''
s='00000000 00000000 00000000 00000000 00000000 10000000 00000000 00000000'
Hight_forward(10,s)