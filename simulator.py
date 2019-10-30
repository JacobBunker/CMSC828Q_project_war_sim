from graphics import *
poly_arr = []
time_arr = []

width = 8
height = 6

objects = 2

def readInput():
    f = open("sim_results.txt")
    data = f.read()
    print(data)
    lines = data.split('\n')
    print(lines)
    i = 0
    while(i < (len(lines) - 1)):
        temp = []
        temp1 = []
        line = lines[i].split()
        print(line)
        time_arr.append(line[0])
        ii = 0
        while(ii < 4):
            temp1.append(Point(line[3+(2*ii)], line[3+(2*ii)+1]))
            ii += 1
        temp.append(temp1)

        temp1 = []
        line = lines[i+1].split()
        print(line)
        ii = 0
        while(ii < 4):
            temp1.append(Point(line[3+(2*ii)], line[3+(2*ii)+1]))
            ii += 1
        temp.append(temp1)
        poly_arr.append(temp)
        i += 2
    i = 0
    for p in poly_arr:
        print( time_arr[i], p[0])
        i += 1

def setWindow():
    win = GraphWin(width=width*100, height=height*100)
    win.setCoords(-25,-100,75,75)
    win.setBackground("green")
    #rect = Rectangle(Point(0.1, 0.1), Point(width-0.1,height-0.1))
    #rect.draw(win)
    #rect.setFill("#654321")
    return win

def setObjects(win, polygons):
    p = Polygon(polygons[1])
    p.draw(win)
    #p.setFill("black")
    p.setOutline("red")

    p = Polygon(polygons[0])
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
        win.getMouse()
        i += 1
    win.close()

#MacOS fix 2
#tk.Toplevel(_root).destroy()

# MacOS fix 1
update()

if __name__ == "__main__":
    simulate()
