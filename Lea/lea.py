from gurobipy import *
from cao_pylib import operation
from cao_pylib import sbox_model


def mand(m,a,b,c,y):
    ta=[]
    tb=[]
    ty=[]
    tc=[]
    for i in range(32):
        ta.append(a[31-i])
        #tb.append([0,0])
        tb.append(b[31-i])
        ty.append(y[31-i])
    for i in range(31):
        tc.append(c[30-i])
    for i in range(32):
        if i==0:
            sbox_model.gen_model_from_ine(m,ta[i]+tb[i]+ty[i],'ine//xor_ine.txt')
            sbox_model.gen_model_from_ine(m,ta[i]+tb[i]+tc[i],'ine//2and_ine.txt')
        elif i<31:
            sbox_model.gen_model_from_ine(m,ta[i]+tb[i]+tc[i-1]+ty[i],'ine//3xor_ine.txt')
            #print(a[i]+b[i]+c[i-1]+y[i])
            sbox_model.gen_model_from_ine(m,ta[i]+tb[i]+tc[i-1]+tc[i],'arx//进位_3_ine.txt')
        else:
            sbox_model.gen_model_from_ine(m,ta[i]+tb[i]+tc[i-1]+ty[i],'ine//3xor_ine.txt')


class Lea():
    def __init__(self,rounds,block_size,index,value):
        self.r=rounds
        self.index=index
        self.block_size=block_size
        if value==0:
            self.vaule=[0,0]
        if value==1:
            self.vaule=[0,1]
        
        self.m=Model()
        self.x=self.m.addVars([i for i in range((self.r+1)*128*2)],vtype=GRB.BINARY,name='x')
        self.c=self.m.addVars([i for i in range(self.r*93*2)],vtype=GRB.BINARY,name='c')

        constr=LinExpr()
        for i in range(256):
            constr+=self.x[i]
        self.m.addConstr(constr>=1)
        for i in range(self.block_size):
            self.m.addConstr(self.x[i*2]+self.x[i*2+1]<=1)
    
    def gen_constr(self):
        for r in range(self.r):
            xindex=r*256
            cindex=r*186

            x0=[] #轮输入
            x1=[] #轮输出
            c=[]  #模加操作的进位
            for i in range(4):
                tmpx0=[]
                tmpx1=[]
                for j in range(32):
                    tmpx0.append([self.x[xindex+(i*32+j)*2],self.x[xindex+(i*32+j)*2+1]])
                    tmpx1.append([self.x[xindex+(i*32+j)*2+256],self.x[xindex+(i*32+j)*2+257]])
                x0.append(tmpx0)
                x1.append(tmpx1)
            for i in range(3):
                tmpc=[]
                for j in range(31):
                    tmpc.append([self.c[cindex+(i*31+j)*2],self.c[cindex+(i*31+j)*2+1]])
                c.append(tmpc)
            self.m.update()
            
            #模加操作的输出
            tmpa=[]
            tmpa.append(operation.roate_right(x1[0],9))
            tmpa.append(operation.roate_left(x1[1],5))
            tmpa.append(operation.roate_left(x1[2],3))

            #模加操作
            for i in range(3):
                mand(self.m,x0[i],x0[i+1],c[i],tmpa[i])
            
            #x0[0]==x1[3]
            for i in range(32):
                for j in range(2):
                    self.m.addConstr(x0[0][i][j]==x1[3][i][j])
            
        index=(self.r*128+self.index)*2
        self.m.addConstr(self.x[index]==self.vaule[0],name='output0')
        self.m.addConstr(self.x[index+1]==self.vaule[1],name='output1')


class Lea_inv():
    def __init__(self,rounds,block_size,index,value):
        self.r=rounds
        self.index=index
        self.block_size=block_size
        if value==0:
            self.value=[0,0]
        if value==1:
            self.value=[0,1]
        
        self.m=Model()
        self.x=self.m.addVars([i for i in range((self.r+1)*256)],vtype=GRB.BINARY,name='x')
        self.c=self.m.addVars([i for i in range(self.r*186)],vtype=GRB.BINARY,name='c')

        constr=LinExpr()
        for i in range(256):
            constr+=self.x[i]
        self.m.addConstr(constr>=1)

        for i in range(self.block_size):
            self.m.addConstr(self.x[i*2]+self.x[i*2+1]<=1)

    def gen_constr(self):
        for r in range(self.r):
            xindex=r*256
            cindex=r*186

            x0=[]
            x1=[]
            c=[]
            for i in range(4):
                tmpx0=[]
                tmpx1=[]
                for j in range(32):
                    tmpx0.append([self.x[xindex+(i*32+j)*2],self.x[xindex+(i*32+j)*2+1]])
                    tmpx1.append([self.x[xindex+(i*32+j)*2+256],self.x[xindex+(i*32+j)*2+257]])
                x0.append(tmpx0)
                x1.append(tmpx1)
            for i in range(3):
                tmpc=[]
                for j in range(31):
                    tmpc.append([self.c[cindex+(i*31+j)*2],self.c[cindex+(i*31+j)*2+1]])
                c.append(tmpc)
            
            #模减操作的输入
            tmps=[]
            tmps.append(operation.roate_right(x0[0],9))
            tmps.append(operation.roate_left(x0[1],5))
            tmps.append(operation.roate_left(x0[2],3))

            #模减操作
            for i in range(3):
                mand(self.m,tmps[i],x1[i],c[i],x1[i+1])
            
            #x0[3]==x1[0]
            for i in range(32):
                for j in range(2):
                    self.m.addConstr(x0[3][i][j]==x1[0][i][j])
        index=(self.r*128+self.index)*2
        
        self.m.addConstr(self.x[index]==self.value[0],name='output0')
        self.m.addConstr(self.x[index+1]==self.value[1],name='output1')
        


def _get_value(lea):
    x=[0 for i in range((lea.r+1)*256)]
    c=[0 for i in range(lea.r*186)]

    for v in lea.m.getVars():
        if v.varName[0]=='x':
            x[eval(v.varName[2:len(v.varName)-1])]=int(v.x)
        if v.varName[0]=='c':
            c[eval(v.varName[2:len(v.varName)-1])]=int(v.x)
    for r in range(lea.r+1):
        print(r)
        s=''
        for i in range(128):
            if i%32==0:
                s+=' '
            if [x[r*256+i*2],x[r*256+i*2+1]]==[0,0]:
                s+='0'
            if [x[r*256+i*2],x[r*256+i*2+1]]==[0,1]:
                s+='1'
            if [x[r*256+i*2],x[r*256+i*2+1]]==[1,1]:
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
        #print()

'''
for i in range(61,62):
    print(i)
    test=Lea_inv(4,i,1)
    test.gen_constr()
    test.m.write('Lea//test1.lp')
    #test.m.Params.OutputFlag=0
    test.m.Params.PoolSolutions=2
    test.m.optimize()

    while test.m.status!=GRB.Status.INFEASIBLE:
        _get_value(test)
        x=[0 for i in range(9999)]
        for v in test.m.getVars():
            if v.varName[0]=='x':
                x[eval(v.varName[2:len(v.varName)-1])]=int(v.x)
        for i in range(128):
            if x[i]==1:
                test.m.addConstr(test.x[i]!=1)
    
'''
test=Lea(3,64,1,1)
test.gen_constr()
test.m.write('test.lp')