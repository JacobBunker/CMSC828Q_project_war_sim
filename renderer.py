from graphics import *
import time

poly_arr = []
time_arr = []

width = 8
height = 6

objects = 5

def readInput():
    f = open("sim_results.txt")
    data = f.read()
    print(data)
    lines = data.split('\n')
    print(lines)
    i = 0
    while(i < (len(lines) - 1)):
        temp = []

        ii = 0
        while(ii < objects):
            temp1 = []
            line = lines[i + ii].split()
            print(line)
            iii = 0
            while(iii < 4):
                temp1.append(Point(line[3+(2*iii)], line[3+(2*iii)+1]))
                iii += 1
            temp.append(temp1)
            ii += 1

        time_arr.append(line[0])

        poly_arr.append(temp)
        i += objects

    i = 0
    for p in poly_arr:
        print( time_arr[i], p[0])
        i += 1

def setWindow():
    win = GraphWin(width=width*100, height=height*100)
    win.setCoords(-500,-500,500,500)
    win.setBackground("green")
    #rect = Rectangle(Point(0.1, 0.1), Point(width-0.1,height-0.1))
    #rect.draw(win)
    #rect.setFill("#654321")
    return win

def setObjects(win, polygons):

    for polygon in polygons:
        p = Polygon(polygon)
        p.draw(win)
        p.setFill("black")
        p.setOutline("red")

def simulate():
    readInput();
    win = setWindow()
    i = 0
    while i < len(poly_arr):
        setObjects(win, poly_arr[i])
        print(time_arr[i])
        time.sleep(0.001)
        i += 1
    win.getMouse()
    win.close()

#MacOS fix 2
#tk.Toplevel(_root).destroy()

# MacOS fix 1
update()

if __name__ == "__main__":
    simulate()
