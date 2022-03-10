from gurobipy import *
from cao_pylib import sbox_model
from cao_pylib import operation

def mand(m,a,b,c,y):
    width=len(a)
    ta=[]
    tb=[]
    ty=[]
    tc=[]
    for i in range(width):
        ta.append(a[width-1-i])
        #tb.append([0,0])
        tb.append(b[width-1-i])
        ty.append(y[width-1-i])
    for i in range(width-1):
        tc.append(c[width-2-i])
    for i in range(width):
        if i==0:
            sbox_model.gen_model_from_ine(m,ta[i]+tb[i]+ty[i],'ine//xor_ine.txt')
            sbox_model.gen_model_from_ine(m,ta[i]+tb[i]+tc[i],'ine//2and_ine.txt')
        elif i<width-1:
            sbox_model.gen_model_from_ine(m,ta[i]+tb[i]+tc[i-1]+ty[i],'ine//3xor_ine.txt')
            #print(a[i]+b[i]+c[i-1]+y[i])
            sbox_model.gen_model_from_ine(m,ta[i]+tb[i]+tc[i-1]+tc[i],'arx//进位_3_ine.txt')
        else:
            sbox_model.gen_model_from_ine(m,ta[i]+tb[i]+tc[i-1]+ty[i],'ine//3xor_ine.txt')


class Speck():
    def __init__(self,rounds,block_size,index,value):
        self.rounds=rounds
        self.index=index
        self.block_size=block_size
        self.half_block=int(self.block_size/2)
        if value==0:
            self.value=[0,0]
        if value==1:
            self.value=[0,1]
        if block_size==32:
            self.shift=[7,2]
        else:
            self.shift=[8,3]
        
        self.m=Model()

        self.x=self.m.addVars([i for i in range(self.block_size*(self.rounds+1)*2)],vtype=GRB.BINARY,name='x')
        self.c=self.m.addVars([i for i in range((self.half_block-1)*self.rounds*2)],vtype=GRB.BINARY,name='c')

        input_diff=LinExpr()
        for i in range(self.block_size*2):
            input_diff+=self.x[i]
        self.m.addConstr(input_diff>=1)

    def gen_constr(self):
        for r in range(self.rounds):
            x0=[]
            x1=[]
            x0_L=[]
            x0_R=[]
            x1_L=[]
            x1_R=[]
            c=[]
            for i in range(self.block_size):
                x0.append([self.x[(r*self.block_size+i)*2],self.x[(r*self.block_size+i)*2+1]])
                x1.append([self.x[((r+1)*self.block_size+i)*2],self.x[((r+1)*self.block_size+i)*2+1]])
            for i in range(self.half_block-1):
                c.append([self.c[(r*(self.half_block-1)+i)*2],self.c[(r*(self.half_block-1)+i)*2+1]])
            
            for i in range(self.half_block):
                x0_L.append(x0[i])
                x0_R.append(x0[i+self.half_block])
                x1_L.append(x1[i])
                x1_R.append(x1[i+self.half_block])

            x0_R_alpha=operation.roate_right(x0_R,self.shift[0])
            x0_L_beta=operation.roate_left(x0_L,self.shift[1])
            
            mand(self.m,x0_R_alpha,x0_R,c,x1_L)
            for i in range(self.half_block):
                
                sbox_model.gen_model_from_ine(self.m,x1_L[i]+x0_L_beta[i]+x1_R[i],'ine//xor_ine.txt')
        
        index=self.index+self.rounds*self.block_size
        self.m.addConstr(self.x[index*2]==self.value[0])
        self.m.addConstr(self.x[index*2+1]==self.value[1])
            

for i in range(64):
    print(i)
    test=Speck(4,64,i,1)
    test.gen_constr()
    test.m.Params.OutputFlag=0
    test.m.optimize()
    if test.m.Status!=GRB.status.INFEASIBLE:
        print(i,'y')