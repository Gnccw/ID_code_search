fp=open('ID.txt','r')
s=fp.read()
fp.close()

s=s.split('\n')

ss=[]
for i in s:
    tmp=''
    
    if len(i)!=0:
        counter=0
        for j in i:
            counter+=1
            if counter<=16:
                if j=='1':
                    tmp+='\?'
                else:
                    tmp+=j
        tmp+=''
        print(tmp)