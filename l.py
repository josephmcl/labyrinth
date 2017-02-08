'''lb.py
    Labyrinth maze generator
'''
from math import acos, pi, ceil, cos, sin, radians
import random
from sys import argv
from PIL import Image, ImageDraw
from time import time

imgx = 4000; imgy = 4000



class Labyrinth:
    def __init__(self, n, d):
        self.n = n                 # number of rings
        self.cuts = d            # number of sections in smallest ring
        self.x = 20 * n            # size of image's X axis
        self.y = 20 * n            # size of image's Y axis
        self.par = 0            # number of nodes created
        self.delta = 360.0/d    # angle of a section between nodes
        self.rep = [self.init_cells(i) for i in range(1,n)]
        # list containing nodes
        self.flat = [item for sublist in self.rep for item in sublist]
        # self.rep flattened
        self.ring_sizes = [len(i) for i in self.rep]
        self.graph = []
        #list containing unsorted edges
        self.order = []
        self.init_edges(self.rep, self.par)
        self.kruskal()

    def init_cells(self, ring):
        '''    arg:        ring:int
            returns     rv:list

            Determine size of ring and return list of unique ints
            with respect to other rings.
        '''
        ang = (180.0*acos((2.0*(ring**2)-1.0)/(2.0*ring**2)))/pi
        if ang <= self.delta / 2:
            self.cuts = self.cuts * 2
            self.delta = ang
        rv = [i+self.par for i in range(self.cuts)]
        self.par = self.par + self.cuts
        return rv

    def cart(self, n, translator):
        '''    args:         n:int, translator:list
            returns:    x:int, y:int

            Takes an int and returns its x,y coord with respect to the
            translator list.
        '''
        r = 1
        for i in range(len(translator)):
            if n >= translator[i]:
                r = r + 1
                n = n - translator[i]
            else:
                break
        r = r*1.0
        t = (n*1.0)/translator[i] * 360
        x = (r+0.5) * cos(radians(t+10))
        y = (r+0.5) * sin(radians(t+10))
        return x*10, y*10

    def find(self, parent, i):
        if parent[i] == i:
            return i
        return self.find(parent, parent[i])

    def union(self, parent, rank, x, y):
        xroot = self.find(parent, x)
        yroot = self.find(parent, y)

        if rank[xroot] < rank[yroot]:
            parent[xroot] = yroot
        elif rank[xroot] > rank[yroot]:
            parent[yroot] = xroot
        else :
            parent[yroot] = xroot
            rank[xroot] += 1

    def kruskal(self):
        result = [] # list containing the final MST

        i = 0 # An index variable, used for sorted edges
        e = 0 # An index variable, used for result[]
        self.graph = sorted(self.graph,key=lambda item: item[2])

        parent =     [i for i in range(self.par)]
        rank =         [0 for i in range(self.par)]

        

        while e < self.par-1 :

            if i >= len(self.graph):
                break
            u,v,w = self.graph[i]
            i = i+1
            x = self.find(parent, u)
            y = self.find(parent ,v)

            if x != y:
                e = e+1
                result.append([u,v,w])
                self.union(parent, rank, x, y)

        for u,v,weight in result:
             self.order.append((u, v))


    def init_edges(self, form, size):
        ''' args:        form:list, size:int
            returns:    -

            Creates possible edges from list of nodes
        '''
        for i in range(len(form)):
            if i < len(form)-1:
                if len(form[i+1]) > len(form[i]): # if next ring is bi-sected
                    for j in range(len(form[i])):
                        self.graph.append([form[i][j], form[i+1][j*2], random.randint(0,20)])
                else:
                    for j in range(len(form[i])):
                        self.graph.append([form[i][j], form[i+1][j], random.randint(0,20)])
            for j in range(len(form[i])):
                self.graph.append([form[i][j], form[i][(j+1) % len(form[i])], random.randint(0,5)])

    def img(self):
        ''' draw!!
        '''
        image = Image.new("RGB", (self.x, self.y),(10,100,10))
        draw = ImageDraw.Draw(image)
        for mv in self.order:

            (a, b) = mv
            x1, y1 = self.cart(a,self.ring_sizes)
            x2, y2 = self.cart(b,self.ring_sizes)
            #print ("({},{}) --> ({},{})").format(x1, y1, x2, y2)
            draw.line((x1+(self.x/2),y1+(self.y/2),x2+(self.x/2),y2+(self.y/2)), fill=(255,255,0), width=4)

        del draw
        image.save('{}.png'.format(argv[1]))
        print("Total nodes:      {}".format(self.par))
        print("Dimensions:       {}x{}".format(self.x, self.y))
        print("Created:          \"{}.png\"".format(argv[1]))



if __name__ == '__main__':
    s = time()
    random.seed("L-A-B-Y-R-I-N-T-H")
    #turtle.speed(speed=0)
    l = Labyrinth(int(argv[2]), int(argv[3]))
    l.img()
    e = time()
    print("Time Elapsed:     [{}]".format(e-s))

