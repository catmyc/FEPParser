class kls:
    def __init__(self, lst):
        self.l = lst

L = list(range(100))

l1 = kls(L)
l2 = kls(L)

del l1.l
print(l2.l)
print(L)

