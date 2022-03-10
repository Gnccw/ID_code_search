from gurobipy import *
from cao_pylib import sbox_model
from cao_pylib import operation
"""
[a127,...,a0  ]    LSB state=256
[a255,...,a128]
[a383,...,a256]
[a511,...,a384]    MSB
"""

def inv_list(x):
    res=[]
    for i in range(len(x)):
        res.append(x[len(x)-1-i])
    return res
class Rectangle():
    def __init__(self,rounds,block_size,index,value):
        self.r=rounds
        self.block_size=block_size
        self.m=Model()
        self.index=index
        if value==0:
            self.value=[0,0]
        if value==1:
            self.value=[0,1]

        self.x=self.m.addVars([i for i in range((self.r+1)*block_size*2)],vtype=GRB.BINARY,name='x')
        self.sb_out=self.m.addVars([i for i in range(self.r*block_size*2)],vtype=GRB.BINARY,name='sb_out')
        self.d=self.m.addVars([i for i in range(self.block_size)],vtype=GRB.BINARY,name='d')

        self.shift_row=[0,1,12,13]
        
        constr=LinExpr()
        for i in range(self.block_size*2):
            constr+=self.x[i]
        self.m.addConstr(constr>=1)

    def gen_constr(self):
        for r in range(self.r):
            #非线性层
            xindex=r*self.block_size*2
            sboxnum=int(self.block_size/4)
            for col in range(sboxnum):
                x=[]
                y=[]
                for row in range(4):
                    index=xindex+2*(row*sboxnum+col)
                    x.append(self.x[index+1])
                    x.append(self.x[index+0])
                    y.append(self.sb_out[index+1])
                    y.append(self.sb_out[index+0])
                
                x=inv_list(x)
                y=inv_list(y)
                self.m.update()
                sbox_model.gen_model_from_ine(self.m,x+y,'rectangle//rectangle_sbox_ine.txt')
            
            #线性层
            for row in range(4):
                x=[]
                y=[]
                for col in range(sboxnum*2):
                    x.append(self.sb_out[xindex+sboxnum*2*row+col])
                    y.append(self.x[xindex+self.block_size*2+sboxnum*2*row+col])
                self.m.update()
                
                sbout=[]
                for i in range(sboxnum):
                    tmp=[0,0]
                    tmp[0]=x[2*i]
                    tmp[1]=x[2*i+1]
                    sbout.append(tmp)
                sbout=operation.roate_left(sbout,self.shift_row[row])
                self.m.update()
                for col in range(sboxnum):
                    self.m.addConstr(sbout[col][0]==y[col*2])
                    self.m.addConstr(sbout[col][1]==y[col*2+1])
        self.index=self.r*self.block_size*2+self.index*2
        self.m.addConstr(self.x[self.index]==self.value[0],name='output0')
        self.m.addConstr(self.x[self.index+1]==self.value[1],name='output1')


class Rectangle_inv():
    def __init__(self,rounds,block_size,index,value):
        self.r=rounds
        self.block_size=block_size
        self.m=Model()
        self.index=index
        if value==0:
            self.value=[0,0]
        if value==1:
            self.value=[0,1]

        self.x=self.m.addVars([i for i in range(self.r*block_size*2)],vtype=GRB.BINARY,name='x')
        self.sb_out=self.m.addVars([i for i in range(self.r*block_size*2)],vtype=GRB.BINARY,name='sb_out')
        self.d=self.m.addVars([i for i in range(self.block_size)],vtype=GRB.BINARY,name='d')

        self.shift_row=[0,1,12,13]
        
        constr=LinExpr()
        for i in range(self.block_size*2):
            constr+=self.x[i]
        self.m.addConstr(constr>=1)

        '''
        #输入差分为2的个数
        constr=LinExpr()
        for i in range(self.block_size):
            constr+=self.d[i]
        self.m.setObjective(constr,GRB.MAXIMIZE)
        for i in range(self.block_size):
            self.m.addConstr(self.x[i*2]>=self.d[i])
            self.m.addConstr(self.x[i*2+1]>=self.d[i])
            self.m.addConstr(self.x[i*2]+self.x[i*2+1]-self.d[i]<=1)
        '''
    
    def gen_constr(self):
        for r in range(self.r-1):
            
            #非线性层
            xindex=r*self.block_size*2
            sboxnum=int(self.block_size/4)
            for col in range(sboxnum):
                x=[]
                y=[]
                for row in range(4):
                    index=xindex+2*(row*sboxnum+col)
                    x.append(self.x[index+1])
                    x.append(self.x[index])
                    y.append(self.sb_out[index+1])
                    y.append(self.sb_out[index])
                x=inv_list(x)
                y=inv_list(y)
                sbox_model.gen_model_from_ine(self.m,x+y,'rectangle//rectangle_sbox_ine.txt')
            
            #线性层
            for row in range(4):
                x=[]
                y=[]
                for col in range(sboxnum*2):
                    x.append(self.sb_out[xindex+sboxnum*2*row+col])
                    y.append(self.x[xindex+self.block_size*2+sboxnum*2*row+col])
                self.m.update()
                
                sbout=[]
                for i in range(sboxnum):
                    tmp=[0,0]
                    tmp[0]=x[2*i]
                    tmp[1]=x[2*i+1]
                    sbout.append(tmp)
                sbout=operation.roate_right(sbout,self.shift_row[row])
                self.m.update()
                for col in range(sboxnum):
                    self.m.addConstr(sbout[col][0]==y[col*2])
                    self.m.addConstr(sbout[col][1]==y[col*2+1])
        #非线性层
        xindex=(self.r-1)*self.block_size*2
        sboxnum=int(self.block_size/4)
        for col in range(sboxnum):
            x=[]
            y=[]
            for row in range(4):
                index=xindex+2*(row*sboxnum+col)
                x.append(self.x[index+1])
                x.append(self.x[index])
                y.append(self.sb_out[index+1])
                y.append(self.sb_out[index])
            x=inv_list(x)
            y=inv_list(y)
            sbox_model.gen_model_from_ine(self.m,x+y,'knot\knot_ine_invert.txt')

        self.index=(self.r-1)*self.block_size*2+self.index*2
        self.m.addConstr(self.sb_out[self.index]==self.value[0],name='output0')
        self.m.addConstr(self.sb_out[self.index+1]==self.value[1],name='output1') 


def get_value(knot,block_size):
    xnum=len(knot.x)
    snum=len(knot.sb_out)
    x=[0 for i in range(xnum)]
    sout=[0 for i in range(snum)]
    rounds=int(snum/(2*block_size))
    for v in knot.m.getVars():
        if v.varName[0]=='x':
            x[eval(v.varName[2:len(v.varName)-1])]=int(v.x)
        if v.varName[0]=='s':
            sout[eval(v.varName[7:len(v.varName)-1])]=int(v.x)
    bx=[]
    bs=[]
    for i in range(int(len(x)/2)):
        tmp=[x[2*i],x[2*i+1]]
        if tmp==[0,0]:
            bx.append(0)
        if tmp==[0,1]:
            bx.append(1)
        if tmp==[1,1]:
            bx.append(2)
    for i in range(int(len(sout)/2)):
        tmp=[sout[2*i],sout[2*i+1]]
        if tmp==[0,0]:
            bs.append(0)
        if tmp==[0,1]:
            bs.append(1)
        if tmp==[1,1]:
            bs.append(2)
    
    xx=''
    ss=''
    for r in range(rounds):
        for row in range(4):
            for col in range(int(block_size/4)):
                xx+=str(bx[r*block_size+row*int(block_size/4)+col])
                ss+=str(bs[r*block_size+row*int(block_size/4)+col])
            xx+='\n'
            ss+='\n'
        xx+='\n'
        ss+='\n'

    if xnum>snum:
        for row in range(4):
            for col in range(int(block_size/4)):
                xx+=str(bx[rounds*block_size+row*int(block_size/4)+col])
            xx+='\n'
        xx+='\n'
    xx=xx.split('\n\n')
    ss=ss.split('\n\n')

    for i in range(rounds):
        print(xx[i])
        print('L')
        print(ss[i])
        print('S')
        print()
    print(xx[rounds])

'''
test=Inv_knot(8,256,185,1)
test.gen_constr()
test.m.optimize()
get_value(test,256)

for i in range(511,512):
    print(i)
    test=Knot(9,512,i,1)
    test.gen_constr()
    test.m.Params.OutputFlag=0
    #test.m.write('test.lp')
    test.m.optimize()
    #print(len(test.sb_out))
    if test.m.status!=GRB.Status.INFEASIBLE:
        print(i,'yes')
        get_value(test,512)
        break
'''  




            

    
    
            

    