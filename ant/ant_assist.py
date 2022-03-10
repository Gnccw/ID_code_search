def perm_128(i):
    if i%4==0:
        return (i+30)%64
    if i%4==1:
        return (i+2)%64
    if i%4==2:
        return (i+58)%64
    if i%4==3:
        return (i+54)%64

def perm_256(i):
    if i%4==0:
        return (i+62)%128
    if i%4==1:
        return (i+2)%128
    if i%4==2:
        return (i+122)%128
    if i%4==3:
        return (i+118)%128

perm128=[]
for i in range(64):
    perm128.append(perm_128(i))
print(perm128)
perm256=[]
for i in range(128):
    perm256.append(perm_256(i))
print(perm256)
