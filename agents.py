from graphics import Circle, Line, Point
import math

theta = math.pi/4
class Fighter():
    def __init__(self, points, win, color):
        self.win = win
        self.radius = 0.1
        self.position = points[0]
        self.looking = points[1]
        self.fire_line = self.get_fire()
        self.color = color
        self.view_left = self.get_left()
        self.view_left.setFill(self.color)
        self.view_right = self.get_right()
        self.view_right.setFill(self.color)
        self.firing = True
        self.view = False
        self.body = Circle(self.position, self.radius)
        self.body.setFill(self.color)


    def draw(self):
        self.body.draw(self.win)
        self.toggle_fire()
        self.toggle_view()
    
    def die(self):
        self.fire_line.undraw()
        self.view_left.undraw()
        self.view_right.undraw()
        self.body.setFill("orange")

    def undraw(self):
        self.body.undraw()

    def update(self, points):
        self.update_position(points[0])
        self.update_view(points[1])

    def update_position(self, position):
        self.body.move(position.x-self.position.x, position.y-self.position.y)
        self.position = position

    def update_view(self, looking):
        self.fire_line.undraw()
        self.view_left.undraw()
        self.view_right.undraw()
        self.looking = looking
        self.fire_line = self.get_fire()
        self.view_left = self.get_left()
        self.view_right = self.get_right()
        if self.firing:
            self.fire_line.draw(self.win)
        if self.looking:
            self.view_left.draw(self.win)
            self.view_right.draw(self.win)


    def toggle_fire(self):
        self.firing = not self.firing
        if self.firing:
            self.fire_line.draw(self.win)
        else:
            self.fire_line.undraw()

    def toggle_view(self):
        self.view = not self.view
        if self.view:
            self.view_left.draw(self.win)
            self.view_right.draw(self.win)
        else:
            self.view_left.undraw()
            self.view_right.undraw()

    def get_left(self):
        x = self.position.x + (self.looking.x - self.position.x)*math.cos(theta) - (self.looking.y - self.position.y)*math.sin(theta)
        y = self.position.y + (self.looking.x - self.position.x)*math.sin(theta) + (self.looking.y - self.position.y)*math.cos(theta)
        l = Line(self.position, Point(x,y))
        l.setFill(self.color)
        return l

    def get_right(self):
        x = self.position.x + (self.looking.x - self.position.x)*math.cos(-theta) - (self.looking.y - self.position.y)*math.sin(-theta)
        y = self.position.y + (self.looking.x - self.position.x)*math.sin(-theta) + (self.looking.y - self.position.y)*math.cos(-theta)
        l = Line(self.position, Point(x,y))
        l.setFill(self.color)
        return l

    def get_fire(self):
        l = Line(self.position, self.looking)
        l.setFill("orange")
        return l


    

    def __repr__(self):
        return "Fighter({}, {})".format(str(self.position), str(self.radius))
