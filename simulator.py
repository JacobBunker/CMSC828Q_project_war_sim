from graphics import *
from agents import *
from sys import stdin
poly_arr = []
#[[Point(2, 3), Point(4, 6), Point(2, 6)], [Point(9, 1), Point(5, 6), Point(3, 8), Point(2, 7)], [Point(6, 2), Point(3.5, 0.2), Point(8.43, 6.73)]]
fight_arr = []
#[[Point(2, 2), Point(3,3), "red"], [Point(8, 8), Point(7,7), "blue"], [Point(10, 2), Point(9,2), "orange"]]
move_list = []
#[[[Point(1.5, 4), Point(1.5, 5)], [Point(6, 8.5), Point(5, 8.5)], [Point(9, 2), Point(8, 2)]]]
            # , 
            # [Point(1.5, 6.25), Point(2, 6.75), Point(4, 8.5), Point(3, 8.5)], 
            # [Point(2.5, 6.25), Point(1, 7.5), Point(2, 8.5), Point(1.5, 7.5)],
            # [Point(2.5, 6.25), Point(1, 7.5), Point(1.5, 7.5), Point(1.5, 6.5)],
            # [Point(2.5, 6.25), Point(1.5, 7.5), Point(1.5, 7.5), Point(2.5, 6.25)]]
spawn_list = []

colors = ["red", "blue", "purple"]

width = 20
height = 20
factor = 30.0
def readInput():
    pos_file = open("positions.txt", "r")
    lines = pos_file.readlines()
    temp = []
    obstacles = True
    for line in lines:
        print(line)
        if line[0] == 'b':
            l = line.split()
            temp.append(Point(float(l[1]), float(l[2])))
        elif line[0] == 's':
            l = line.split()
            spawn_list.append([Point(float(l[1]),float(l[2])),float(l[3])])
        elif line[0] == 'e':
            poly_arr.append(temp)
            temp =[]
        elif obstacles:
            poly_arr.append(temp)
            l = line.split()
            point = Point(float(l[1]), float(l[2]))
            angle = float(l[2])
            temp = []
            temp.append([point, angle])
            print(poly_arr)
            obstacles = False
        elif line[0] == '0':
            move_list.append(temp)
            temp = []
            print(move_list)
            l = line.split()
            point = Point(float(l[1]), float(l[2]))
            angle = float(l[2])
            temp.append([point, angle])
        else:
            l = line.split()
            point = Point(float(l[1]), float(l[2]))
            angle = float(l[2])
            temp.append([point, angle])
    i = 0
    for start in move_list[0]:
        fight_arr.append([start[0], start[1], colors[i]])
        i += 1


def setWindow():
    win = GraphWin(width=width*factor, height=height*factor)
    win.setCoords(0,0,width,height)
    win.setBackground("#654321")
    return win



def setObstacles(win, polygons):
    for polygon in polygons:
        p = Polygon(polygon)
        p.draw(win)
        p.setFill("green")
        p.setOutline("green")

def setSpawns(win):
    for spawn in spawn_list:
        c = Circle(spawn[0], spawn[1])
        c.draw(win)
        c.setFill("purple")
        c.setOutline("purple")

def setFighters(win, fs):
    fighters = []
    for f in fs:
        fighter = Fighter(f[0:2], win, f[2])
        fighter.draw()
        fighters.append(fighter)
    return fighters

def end_game(win, fighters, winner):
    victor = ""
    for fighter in fighters:
        if fighter is not winner:
            fighter.die()
        else:
            victor = fighter.color
    win.getMouse()
    for fighter in fighters:
        if fighter is not winner:
            fighter.undraw()
        else:
            fighter.toggle_fire()
    t = Text(Point(width/2,height/2), victor +" Fighter Wins")
    t.setSize(30)
    t.draw(win)
    win.getMouse()

def simulate():
    readInput()
    win = setWindow()
    setObstacles(win, poly_arr)
    setSpawns(win)


    fighters = setFighters(win, fight_arr)
    winner = fighters[0]
   
    win.getMouse()
    while len(move_list) > 0:
        moves = move_list.pop(0)
        for fighter in fighters:
            fighter.update(moves[0])
            moves.pop(0)
        # if len(move_list) == 0:
        #     f1.toggle_fire()
        #     f1.toggle_view()
        win.getMouse()
    end_game(win, fighters, winner)
    win.close()

#MacOS fix 2
#tk.Toplevel(_root).destroy()

# MacOS fix 1
update()

if __name__ == "__main__":
    simulate()