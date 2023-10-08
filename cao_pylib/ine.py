import copy

'''
check that whether the point [x_0,x_1,...,x_{n-1}] satisfies the inequality "a_0 * x_0 + a_1 * x_1 + ... + a_{n-1} * x_{n-1} + c >= 0"
'''
def is_sat(point, ineq):
    temp = 0
    for i in range(len(point)):
        temp = temp + point[i] * ineq[i]
    temp = temp + ineq[len(point)]
    if temp < 0:
        return 0
    elif temp >= 0:
        return 1

# return points that don't satisfy the inequality
def collect_cutoff(point, ineq):
    temp = []
    for p in point:
        flag = is_sat(p,ineq)
        if flag == 0:
            temp.append(p)
    return temp

# select N inequalities from the inequalities generated from sage to cutoff impossible points
def greedySelection(N, ine_convex_hull, im_point):
    ine = copy.deepcopy(ine_convex_hull)
    point = copy.deepcopy(im_point)
    select_ine = []
    for i in range(N):
        cutoff = []
        count_of_cutoff = []
        for l in ine:
            cutoff_of_ine = collect_cutoff(point, l)
            cutoff.append(cutoff_of_ine)
            count_of_cutoff.append(len(cutoff_of_ine))
        max_count_index = count_of_cutoff.index(max(count_of_cutoff))
        select_ine.append(ine[max_count_index])
        ine.remove(ine[max_count_index])
        for p in cutoff[max_count_index]:
            point.remove(p)
    return select_ine

# return the inequalities that cutoff all impossible points
def	best_cutoff_ine(ine_convex_hull,im_point):
    ine = copy.deepcopy(ine_convex_hull)
    point = copy.deepcopy(im_point)
    select_ine = []
    while point != []:
        cutoff = []
        count_of_cutoff = []
        for l in ine:
            cutoff_of_ine = collect_cutoff(point, l)
            cutoff.append(cutoff_of_ine)
            count_of_cutoff.append(len(cutoff_of_ine))
        max_count_index = count_of_cutoff.index(max(count_of_cutoff))
        select_ine.append(ine[max_count_index])
        ine.remove(ine[max_count_index])
        for p in cutoff[max_count_index]:
            point.remove(p)
    return select_ine


def list2ine(a):
    """
    transform a list to a ine
    """
    s='('
    for i in range(len(a)-3):
        s+=str(a[i])+', '
    s+=str(a[len(a)-3])+') x '
    if a[len(a)-2] >=0:
        s+='+ '+str(a[len(a)-2])
    else:
        s+=' - '+str(a[len(a)-2]*(-1))
    s+=' >= '+str(a[len(a)-1])
    return s


def clause2ine(a):
    """
    transform a clause to a ine,clause form[0,1,-1,0]->[b+c']
    """
    ine=copy.deepcopy(a)
    const=-1
    for i in range(len(a)):
        if a[i]==-1:
            const+=1
    ine+=[const,0]
    return list2ine(ine)

