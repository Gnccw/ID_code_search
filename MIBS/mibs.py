from gurobipy import *
from cao_pylib import sbox_model

def _xor(m,x,y):
    varnum=len(x)
    for i in range(4):
        tmp=[]
        for j in range(varnum):
            tmp+=x[j][i]
        tmp+=y[i]
        sbox_model.gen_model_from_ine(m,tmp,'ine//_'+str(varnum)+'xor_ine.txt')

def list2bin(x):
    if x==[0,0]:
        return '0'
    if x==[0,1]:
        return '1'
    if x==[1,1]:
        return '2'
        
def _get_value(Mibs):
    xnum=len(Mibs.x)
    rounds=int(xnum/128)

    x=[0 for i in range(xnum)]
    sout=[0 for i in range((rounds-1)*64)]
    xorout=[0 for i in range((rounds-1)*64)]
    for v in Mibs.m.getVars():
        if v.varName[0]=='x' and v.varName[1]=='[':
            x[eval(v.varName[2:len(v.varName)-1])]=int(v.x)
        if v.varName[0]=='s':
            sout[eval(v.varName[5:len(v.varName)-1])]=int(v.x)
        if v.varName[0]=='x' and v.varName[1]=='o':
            xorout[eval(v.varName[7:len(v.varName)-1])]=int(v.x)

           
    for i in range(rounds):
        index_start=i*128
        xx=''
        ssout=''
        xxout=''
        for j in range(64):
            tmp=[x[index_start+j*2],x[index_start+j*2+1]]
            if tmp==[0,0]:
                xx+='0'
            if tmp==[0,1]:
                xx+='1'
            if tmp==[1,1]:
                xx+='2'
        if i<rounds-1:
            for j in range(32):
                tmp=[sout[i*64+j*2],sout[i*64+j*2+1]]
                tmpx=[xorout[i*64+j*2],xorout[i*64+j*2+1]]
                ssout+=list2bin(tmp)
                xxout+=list2bin(tmpx)
        print(xx)
        print(ssout)
        print(xxout)
        print()

class Mibs():
    def __init__(self,rounds,index,value,inv):
        self.r=rounds
        self.index=index
        #inv=1 逆向
        self.inv=inv

        if value==1:
            self.value=[0,1]
        if value==0:
            self.value=[0,0]
        
        self.m=Model()
        self.x=self.m.addVars([i for i in range((self.r+1)*64*2)],vtype=GRB.BINARY,name='x')
        self.sout=self.m.addVars([i for i in range(self.r*32*2)],vtype=GRB.BINARY,name='sout')
        self.xorout=self.m.addVars([i for i in range(self.r*32*2)],vtype=GRB.BINARY,name='xorout')

        self.xorpos=[[1, 2, 3, 4, 5, 6],
                     [0, 2, 3, 5, 6, 7],
                     [0, 1, 3, 4, 6, 7],
                     [0, 1, 2, 4, 5, 7],
                     [0, 1, 3, 4, 5],   
                     [0, 1, 2, 5, 6],   
                     [1, 2, 3, 6, 7],   
                     [0, 2, 3, 4, 7]]

        constr=LinExpr()
        for i in range(128):
            constr+=self.x[i]
        self.m.addConstr(constr>=1)
    
    def gen_constr(self):
        for r in range(self.r):
            x_index_start=r*64*2
            s_index_start=r*32*2

            si=[]
            so=[]
            xo=[]
            for i in range(8):
                tmpsi=[]
                tmpso=[]
                tmpxo=[]
                for j in range(4):
                    tmpsi.append([self.x[x_index_start+(i*4+j)*2],self.x[x_index_start+(i*4+j)*2+1]])
                    tmpso.append([self.sout[s_index_start+(i*4+j)*2],self.sout[s_index_start+(i*4+j)*2+1]])
                    tmpxo.append([self.xorout[s_index_start+(i*4+j)*2],self.xorout[s_index_start+(i*4+j)*2+1]])
                si.append(tmpsi)
                so.append(tmpso)
                xo.append(tmpxo)
            self.m.update()
            #非线性层
            for i in range(8):
                tmp=[]
                for j in range(4):
                    tmp+=si[i][j]
                for j in range(4):
                    tmp+=so[i][j]
                sbox_model.gen_model_from_ine(self.m,tmp,'MIBS//mibs_sbox_ine.txt')
            
            #线性层
            for i in range(8):
                x=[]
                for j in range(len(self.xorpos[i])):
                    x.append(so[self.xorpos[i][j]])
                _xor(self.m,x,xo[i])
            
            x=[]
            y=[]
            z=[]
            for i in range(32):
                x.append([self.xorout[s_index_start+i*2],self.xorout[s_index_start+i*2+1]])
                y.append([self.x[x_index_start+(i+32)*2],self.x[x_index_start+(i+32)*2+1]])
                z.append([self.x[x_index_start+(64+i)*2],self.x[x_index_start+(64+i)*2+1]])
            for i in range(32):
                sbox_model.gen_model_from_ine(self.m,x[i]+y[i]+z[i],'ine//xor_ine.txt')
                #print(x[i]+y[i]+z[i])
            for i in range(64):
                self.m.addConstr(self.x[x_index_start+i]==self.x[x_index_start+i+192])
        index=self.r*128+self.index*2
        self.m.addConstr(self.x[index]==self.value[0])
        self.m.addConstr(self.x[index+1]==self.value[1])
        
for i in range(64):
    print(i)
    test=Mibs(3,i,1,0)
    test.gen_constr()
    test.m.write('MIBS//test.lp')
    
    test.m.Params.OutputFlag=0
    test.m.optimize()
    if test.m.status!=GRB.Status.INFEASIBLE:
        print(i,'yes')
        _get_value(test)
        break 
     