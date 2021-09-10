import matplotlib.pyplot as plt
import random
	
numberOfNodes = 100

k = 1.2
d = 14

r = k * d

points = []
xvalues = []
yvalues = []

def getDistanceOfPointsFromPoint(point):
	distances = []
	for otherPoint in points:
		x1 = point[0]
		y1 = point[1]
		x2 = otherPoint[0]
		y2 = otherPoint[1]
		distance = ((((x2 - x1 )**2) + ( (y2-y1)**2) )**0.5)
		distances.append(distance)
	return distances

def drawCircles():
	for point in points:
		circle = plt.Circle( (point[0], point[1]), r, fill = False)
		plt.gca().add_artist(circle)

def getNeighborNodes(point):
	neighbors = []
	distances = getDistanceOfPointsFromPoint(point)
	for distance in distances:
		if (distance < r and( distance > 0)):
			neighboring = points[distances.index(distance)]
			neighbors.append(neighboring)
	return neighbors

def printTable():
    table = []
    print ("Sensor Index |  Neighbor Indices")
    for point in points:
        neighbors = getNeighborNodes(point)
        indices = []
        for neighbor in neighbors:
            indices.append(points.index(neighbor))
        table.append([points.index(point), indices])
        print(str(points.index(point)) + " | " + str(indices))
	#print(tabulate(table, headers=tableHeaders, tablefmt="grid"))

def generateNodes():
	for i in range(numberOfNodes):
		x = random.randint(0,150)
		y = random.randint(0,150)
		point = [x, y]
		points.append(point)
		xvalues.append(x)
		yvalues.append(y)

def connectNeighbors():
    for point in points:
        neighbors = getNeighborNodes(point)
        for neighbor in neighbors:
            x1 = point[0]
            y1 = point[1]
            x2 = neighbor[0]
            y2 = neighbor[1]
            plt.plot([x1,x2], [y1,y2], color='b')
		
generateNodes()

plt.scatter(xvalues, yvalues)

plt.axis([0,200,0,200])

drawCircles()

printTable()

connectNeighbors()

plt.xlabel('x-axis')

plt.ylabel('y-axis')

plt.show()