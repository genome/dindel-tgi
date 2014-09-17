import math

def addLogs(l1, l2):
    if l1>l2:
        diff = l2 - l1;
        return l1+math.log(1.0+math.exp(diff))
    else:
        diff = l1 - l2;
        return l2 + math.log(1.0+math.exp(diff))

