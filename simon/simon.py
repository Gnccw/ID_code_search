from gurobipy import *
from cao_pylib import operation
from cao_pylib import sbox_model
import time
rounds=18
block_size=128 
class Simon():
    def __init__(self,rounds,blocksize,index,value):
        self.r=rounds
        self.block_size=blocksize
        self.index=index
        if value==1:
            self.value=[0,1]
        if value==0:
            self.value=[0,0]
        
        self.m=Model()
        self.x=self.m.addVars([i for i in range(self.block_size*(self.r+1)*2)],vtype=GRB.BINARY,name='x')
        self.aout=self.m.addVars([i for i in range(self.r*self.block_size)],vtype=GRB.BINARY,name='aout')

        
        #输入差分不为0
        constr=LinExpr()
        for i in range(self.block_size*2):
            constr+=self.x[i]
        self.m.addConstr(constr>=1)

        #限制输入差分为0,1
        for i in range(self.block_size):
            self.m.addConstr(self.x[i*2]==0)

        '''
        #右边为0
        if self.inv==0:
            for i in range(self.block_size):
                self.m.addConstr(self.x[i+blocksize]==0)

        
        if self.inv==1:
        #左边为0
            for i in range(self.block_size):
                self.m.addConstr(self.x[i]==0)
        '''

    def gen_constr(self):
        for r in range(self.r):
            index_start=r*self.block_size*2
            aindex_start=r*self.block_size
            x0=[]
            x1=[]
            x2=[]
            x3=[]
            x4=[]
            x5=[]
            a=[]
            for i in range(int(self.block_size/2)):
                tmp=[self.x[index_start+2*i],self.x[index_start+2*i+1]]
                tmp1=[self.x[index_start+self.block_size+2*i],self.x[index_start+self.block_size+2*i+1]]
                tmp2=[self.x[index_start+self.block_size*2+2*i],self.x[index_start+self.block_size*2+2*i+1]]
                tmp3=[self.aout[aindex_start+i*2],self.aout[aindex_start+i*2+1]]
                x0.append(tmp)
                x1.append(tmp1)
                x5.append(tmp2)
                a.append(tmp3)
            x2=operation.roate_left(x0,1)
            x3=operation.roate_left(x0,8)
            x4=operation.roate_left(x0,2)
            self.m.update()
            
            for i in range(int(self.block_size/2)):
                sbox_model.gen_model_from_ine(self.m,x2[i]+x3[i]+a[i],'ine//2and_ine.txt')
                sbox_model.gen_model_from_ine(self.m,a[i]+x1[i]+x4[i]+x5[i],'ine//3xor_ine.txt')
            
            for i in range(int(self.block_size/2)):
                self.m.addConstr(self.x[index_start+i*2]==self.x[index_start+3*self.block_size+2*i])
                self.m.addConstr(self.x[index_start+i*2+1]==self.x[index_start+3*self.block_size+2*i+1])
    
        
        self.m.addConstr(self.x[self.r*self.block_size*2+self.index*2]==self.value[0],name='output0')
        self.m.addConstr(self.x[self.r*self.block_size*2+self.index*2+1]==self.value[1],name='output1')
        
class Simon_inv():
    def __init__(self,rounds,blocksize,index,value):
        self.r=rounds
        self.block_size=blocksize
        self.index=index
        if value==1:
            self.value=[0,1]
        if value==0:
            self.value=[0,0]
        
        self.m=Model()
        self.x=self.m.addVars([i for i in range(self.block_size*(self.r+1)*2)],vtype=GRB.BINARY,name='x')
        self.aout=self.m.addVars([i for i in range(self.r*self.block_size)],vtype=GRB.BINARY,name='aout')
        self.d=self.m.addVars([i for i in range(self.block_size)],vtype=GRB.BINARY,name='d')

        '''
        #限制输入差分为0,1
        for i in range(self.block_size):
            self.m.addConstr(self.x[i*2]==0)
        '''
        #输入差分不为0
        constr=LinExpr()
        for i in range(self.block_size*2):
            constr+=self.x[i]
        self.m.addConstr(constr>=1)
        
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
        if self.inv==1:
        #左边为0
            for i in range(self.block_size):
                self.m.addConstr(self.x[i]==0)
        '''

    def gen_constr(self):
        for r in range(self.r):
            index_start=r*self.block_size*2
            aindex_start=r*self.block_size
            x0=[]
            x1=[]
            x2=[]
            x3=[]
            x4=[]
            x5=[]
            a=[]
            for i in range(int(self.block_size/2)):
                tmp=[self.x[index_start+2*i],self.x[index_start+2*i+1]]
                tmp1=[self.x[index_start+self.block_size+2*i],self.x[index_start+self.block_size+2*i+1]]
                tmp2=[self.x[index_start+self.block_size*2+2*i],self.x[index_start+self.block_size*2+2*i+1]]
                tmp3=[self.aout[aindex_start+i*2],self.aout[aindex_start+i*2+1]]
                x0.append(tmp)
                x1.append(tmp1)
                x5.append(tmp2)
                a.append(tmp3)
            x2=operation.roate_left(x0,1)
            x3=operation.roate_left(x0,8)
            x4=operation.roate_left(x0,2)
            self.m.update()
            
            for i in range(int(self.block_size/2)):
                sbox_model.gen_model_from_ine(self.m,x2[i]+x3[i]+a[i],'ine//2and_ine.txt')
                sbox_model.gen_model_from_ine(self.m,a[i]+x1[i]+x4[i]+x5[i],'ine//3xor_ine.txt')
            
            for i in range(int(self.block_size/2)):
                self.m.addConstr(self.x[index_start+i*2]==self.x[index_start+3*self.block_size+2*i])
                self.m.addConstr(self.x[index_start+i*2+1]==self.x[index_start+3*self.block_size+2*i+1])
        
        
        self.index=(self.index+int(self.block_size/2))%self.block_size
        self.m.addConstr(self.x[self.r*self.block_size*2+self.index*2]==self.value[0],name='output0')
        self.m.addConstr(self.x[self.r*self.block_size*2+self.index*2+1]==self.value[1],name='output1')


def get_value(simon,block_size):
    xnum=len(simon.x)
    anum=len(simon.aout)
    x=[0 for i in range(xnum)]
    aout=[0 for i in range(anum)]
    for v in simon.m.getVars():
        if v.varName[0]=='x':
            x[eval(v.varName[2:len(v.varName)-1])]=int(v.x)
        if v.varName[0]=='a':
            aout[eval(v.varName[5:len(v.varName)-1])]=int(v.x)
    s=''
    for r in range(simon.r+1):
        for i in range(block_size):
            if [x[(r*simon.block_size+i)*2],x[(r*simon.block_size+i)*2+1]]==[0,0]:
                s+='0'
            if [x[(r*simon.block_size+i)*2],x[(r*simon.block_size+i)*2+1]]==[0,1]:
                s+='1'
            if [x[(r*simon.block_size+i)*2],x[(r*simon.block_size+i)*2+1]]==[1,1]:
                s+='2'
        s+='\n'
    print(s)


def get_inputdiff(cla):
    res=[]
    x=[0 for i in range(len(cla.x))]
    for v in cla.m.getVars():
        if v.varName[0]=='x':
            x[eval(v.varName[2:len(v.varName)-1])]=int(v.x)
  
    for i in range(cla.block_size*2):
        res.append(x[i])
    return res

def get_output_diff(cla):
    res=[]
    x=[0 for i in range(len(cla.x))]
    for v in cla.m.getVars():
        if v.varName[0]=='x':
            x[eval(v.varName[2:len(v.varName)-1])]=int(v.x)
    for i in range(cla.block_size*2):
        res.append(x[i+cla.r*cla.block_size*2])
    return res


#给输入添加约束
def add_input_constr(cla,inputdiff):
    constr=LinExpr()
    const=1
    for i in range(len(inputdiff)):
        if inputdiff[i]==1:
            constr-=cla.x[i]
            const-=1
        else:
            constr+=cla.x[i]
    cla.m.addConstr(constr>=const)

def add_input_value(cla,inputdiff):
    for i in range(len(inputdiff)):
        cla.m.addConstr(cla.x[i]==inputdiff[i])


def add_output_constr(cla,index,s):
    index_1=(index+int(cla.block_size/2))%cla.block_size
    value=[]
    if s==0:
        value=[0,0]
    if s==1:
        value=[0,1]
    cla.m.addConstr(cla.x[cla.r*cla.block_size*2+index_1*2]==value[0],name='output0')
    cla.m.addConstr(cla.x[cla.r*cla.block_size*2+index_1*2+1]==value[1],name='output1')


#打印list
def list2str(inputlist):
    s=''
    for i in range(int(len(inputlist)/2)):
        if [inputlist[i*2],inputlist[i*2+1]]==[0,0]:
            s+='0'
        if [inputlist[i*2],inputlist[i*2+1]]==[0,1]:
            s+='1'
        if [inputlist[i*2],inputlist[i*2+1]]==[1,1]:
            s+='2'
    return s


        

'''
aa=[]
for i in range(64,128):
    print(i)
    snum=0
    inputs=[]
    test=Simon(7,128,i,1,1)
    test.gen_constr()
    #test.m.write('48.lp')
    test.m.Params.OutputFlag=0
    test.m.optimize()
    while test.m.status!=GRB.Status.INFEASIBLE:
        get_value(test,128)
        list2str(get_inputdiff(test))
        add_input_constr(test,snum,get_inputdiff(test))
        test.m.Params.OutputFlag=0
        test.m.optimize()
        snum+=1
        print(snum)
    aa.append(snum)
    fp=open('simon128.txt','a')
    fp.write(str(snum)+'\n')
    fp.close()
print(aa)
#test.m.write('simon32.lp')
'''

#返回所有R轮的矛盾
def  all_poss_contradict(R,block_size):
    res=[]
    R0_R1={32:[5,6],48:[6,7],64:[7,8],96:[9,10],128:[11,12]}
    R0=R0_R1[block_size][0]
    R1=R0_R1[block_size][1]
    print(R0,R1)
    for r in range(R-R1,R0+1):
        print(r)
        for index in range(block_size):
            for value in range(2):
                tmp=[]
                M0=Simon(r,block_size,index,value,0)
                M1=Simon(R-r,block_size,index,value^1,1)
                M0.gen_constr()
                M1.gen_constr()
                M0.m.Params.OutputFlag=0
                M1.m.Params.OutputFlag=0
                M0.m.optimize()
                M1.m.optimize()
                if M0.m.Status!=GRB.status.INFEASIBLE and M1.m.Status!=GRB.status.INFEASIBLE:
                    tmp=[r,index,value]
                    res.append(tmp)
    s=''
    for var in res:
        s+=str(var[0])+','
        s+=str(var[1])+','
        s+=str(var[2])+'\n'
    fp=open(str(block_size)+'_r'+str(R)+'_constra.txt','w')
    fp.write(s)
    fp.close()

'''
all_poss_contradict(10,32)
all_poss_contradict(11,48)
all_poss_contradict(12,64)
all_poss_contradict(15,96)
all_poss_contradict(18,128)
'''
def is_output_diff(diff):
    res=0
    for i in diff:
        if i=='1':
            res+=1
    if res!=1:
        return 0
    else:
        return 1


def all_imposs_diff(rounds,block_size):

    #-----------read constradiction from file------------------
    constra=[]
    fp=open(str(block_size)+'_r'+str(rounds)+'_constra.txt','r')
    s=fp.read()
    fp.close()
    s=s.split('\n')[:len(s.split('\n'))-1]

    for v in s:
        v1=v.split(',')
        tmp=[]
        for i in range(3):
            tmp.append(eval(v1[i]))
        constra.append(tmp)
    #-----------read constradiction from file------------------
    rounds_M0=[100,0]
    for i in range(len(constra)):
        if constra[i][0]<rounds_M0[0]:
            rounds_M0[0]=constra[i][0]
        if constra[i][0]>rounds_M0[1]:
            rounds_M0[1]=constra[i][0]

    all_inputs=[]
    input_output_num=[]

    for v in constra:
        if v[0]!=99:
            print(v)

            #--------------construct M0 and add all_inputs constraint-------------
            M0=Simon(v[0],block_size,v[1],v[2],0)
            if len(all_inputs)!=0:
                for inputs in all_inputs:
                    add_input_constr(M0,inputs)
            M0.gen_constr()
            M0.m.Params.OutputFlag=0
            M0.m.optimize()
            #--------------construct M0 and add all_inputs constraint-------------

            #当M0轮数较小时,输入不可能差分会很多，所以先排除掉不存在不可能输出差分的情形
            flag=0
            M2=Simon(rounds-v[0],block_size,v[1],v[2]^1,1)
            M2.gen_constr()
            M2.m.Params.OutputFlag=0
            M2.m.optimize()
            if M2.m.Status!=GRB.INFEASIBLE:
                flag=1
            #------search all input differential of the given constradiction-----
            while M0.m.Status!=GRB.status.INFEASIBLE and flag==1:
                s='inputdiff:'
                inputdiff_M0=get_inputdiff(M0)
                #get_value(M0,block_size)
                all_inputs.append(inputdiff_M0)
            
                s+=list2str(inputdiff_M0)+'\n'
                add_input_constr(M0,inputdiff_M0)

                M0.m.update()
                #M0.m.write('test.lp')
                M0.m.optimize()

                all_outputdiff_given_inputdiff=[]
                for r in range(rounds_M0[0],rounds_M0[1]+1):
                    M0_1=Simon(r,block_size,1,0,0)
                    M0_1.gen_constr()
                    M0_1.m.update()
                    M0_1.m.remove(M0_1.m.getConstrByName('output0'))
                    M0_1.m.remove(M0_1.m.getConstrByName('output1'))
                    add_input_value(M0_1,inputdiff_M0)
                    M0_1.m.update()
                    M0_1.m.Params.OutputFlag=0
                    M0_1.m.optimize()
                    #get_value(M0_1,block_size)
                    

                    outputdiff_M0_1=get_output_diff(M0_1)
                    outputdiff_str_M0_1=list2str(outputdiff_M0_1)

                    for i in range(block_size):
                        if outputdiff_str_M0_1[i]!='2':
                            if [r,i,eval(outputdiff_str_M0_1[i])] in constra:
                                M1=Simon(rounds-r,block_size,i,eval(outputdiff_str_M0_1[i])^1,1)
                                M1.gen_constr()

                                if len(all_outputdiff_given_inputdiff)!=0:
                                    for i in range(len(all_outputdiff_given_inputdiff)):
                                        add_input_constr(M1,all_outputdiff_given_inputdiff[i])
                                
                                M1.m.Params.OutputFlag=0
                                M1.m.optimize()

                                while M1.m.Status!=GRB.status.INFEASIBLE:
                                    #get_value(M1,block_size)
                                    imposs_output_diff=get_inputdiff(M1)
                                    s+='          '+list2str(imposs_output_diff)+'\n'
                                    all_outputdiff_given_inputdiff.append(imposs_output_diff)
                                    
                                    add_input_constr(M1,imposs_output_diff)
                                    M1.m.update()
                                    M1.m.optimize()
                input_output_num.append(len(all_outputdiff_given_inputdiff))
                s+='\n'
                fp=open(str(block_size)+'_all_imposs_diff_'+str(rounds)+'.txt','a')
                fp.write(s)
                fp.close()
    numss=''
    for i in range(len(input_output_num)):
        numss+=str(input_output_num[i])+','
    fp=open(str(block_size)+'_all_imposs_diff_'+str(rounds)+'.txt','a')
    fp.write(numss[:len(numss)-1]+'\n')
    fp.close()
            

'''
t1=time.time()
all_imposs_diff(rounds,block_size)
t2=time.time()
t=(int((t2-t1)*100))/100
fp=open(str(block_size)+'_all_imposs_diff_'+str(rounds)+'.txt','a')
fp.write('time:'+str(t)+'s')
fp.close()
print(t)




def get_max_output(constra,block_size):
    M0=Simon(constra[0],blocksize,constra[1],constra[2],0)
    for i in range()
'''
'''
blocksize=32
for i in [1,blocksize-3]:
    print(i)
    test=Simon(6,blocksize,i,0,0) #(rounds,blocksize,index,value,inv)
    test.gen_constr()
    test.m.Params.OutputFlag=0
    test.m.optimize()
    if test.m.status!=GRB.Status.INFEASIBLE:
        get_value(test,blocksize)
        #print('1111111')
        #break
'''
'''
test=Simon(7,128,3,0,0)
test.gen_constr()
test.m.optimize()
get_value(test,128)
'''
test=Simon(2,32,1,1)
test.gen_constr()
test.m.optimize()
get_value(test,32)