import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb

defaultSquareColor = [130/360, 0, 100/100]
targetSquareColor = [30/360, 100/100, 100/100]
currentSquareColor = [340/360, 95/100, 100/100]
traversedSquareColor = [200/360, 95/100, 100/100]

global currentSquare
currentSquare = (0,0)

UP = 0
DOWN = 1
RIGHT = 2
LEFT = 3
actions = ['UP', 'DOWN', 'RIGHT', 'LEFT']

numStates = 5 * 5
numActions = 4

curEpisode = 0

meanRewards = []

traversedSquares = []


grid = defaultSquareColor * np.ones((5, 5, 3))
grid[4, 4] = targetSquareColor # color objective square

fig, ax = plt.subplots(figsize=(12, 4))
ax.grid(which='minor')
ax.set_xticks(np.arange(5))
ax.set_xticks(np.arange(5) - 0.5, minor=True)
ax.set_yticks(np.arange(5))
ax.set_yticks(np.arange(5) - 0.5, minor=True)

def getPositionFromId(id):
        return (id // 5), (id % 5)

qText =  [ax.text(*getPositionFromId(i)[::-1], '0', fontsize=11, verticalalignment='center',  horizontalalignment='center') for i in range(5 * 5)]

im = ax.imshow(hsv_to_rgb(grid), interpolation='nearest', vmin=0, vmax=1) 

qValues = np.zeros((numStates, numActions))

df = pd.DataFrame(qValues, columns=[' up ', 'down', 'right', 'left'])

df.index.name = 'States'

df.head()

def reset():
    global currentSquare
    currentSquare = (0, 0)  
    global curEpisode 
    curEpisode += 1      
    global traversedSquares
    traversedSquares.clear()
    traversedSquares.append((0,0))
    return getIdFromPosition(currentSquare)
    
def preformAction(action):
    # Checking for valid options for movement
    global currentSquare, grid
    if action == 0 and currentSquare[0] > 0:
        currentSquare = (currentSquare[0] - 1, currentSquare[1])
    if action == 1 and currentSquare[0] < 4:
        currentSquare = (currentSquare[0] + 1, currentSquare[1])
    if action == 2 and currentSquare[1] < 4:
        currentSquare = (currentSquare[0], currentSquare[1] + 1)
    if action == 3 and currentSquare[1] > 0:
        currentSquare = (currentSquare[0], currentSquare[1] - 1)
            
    # give reward based on current square
    if all(grid[currentSquare] == targetSquareColor):
        reward = 100
        done = True
    else:
        reward = -1
        done = False
    return getIdFromPosition(currentSquare), reward, done

def getIdFromPosition(pos):
        return pos[0] * 5 + pos[1]
    
def drawToGrid(qValues=None, action=None, qMax=False, shouldColor=False):
       
        if shouldColor:      
            grid = defaultSquareColor * np.ones((5, 5, 3))
            grid[4, 4] = targetSquareColor # color objective square
        else:            
            grid = grid.copy()
            
        for traversedSquare in traversedSquares:
            grid[traversedSquare] = traversedSquareColor

        grid[currentSquare] = currentSquareColor

        traversedSquares.append(currentSquare)

        im.set_data(hsv_to_rgb(grid))
               
        if qValues is not None:
            xs = np.repeat(np.arange(5), 5)
            ys = np.tile(np.arange(5), 5)  
            
            for i, text in enumerate(qText):
                if qMax:
                    q = max(qValues[i])    
                    txt = '{:.2f}'.format(q)
                    text.set_text(txt)
                else:                
                    actions = ['U', 'D', 'R', 'L']
                    txt = '\n'.join(['{}: {:.2f}'.format(k, q) for k, q in zip(actions, qValues[i])])
                    text.set_text(txt)
                
        if action is not None:
            global curEpisode
            ax.set_title("Episode: " + str(curEpisode), color='r', weight='bold', fontsize=32)

        plt.pause(0.01)


def greedPolicy(qValues, state, explorationRate):
    if np.random.random() < explorationRate:
        return np.random.choice(4)
    else:
        return np.argmax(qValues[state])

def qLearning(grid, nEpisodes=30, shouldDrawToGrid=True, explorationRate=0.1, learnRate=0.5, gamma=0.9):    
    qValues = np.zeros((numStates, numActions))
    episodeRewards = []
    
    for _ in range(nEpisodes):
        state = reset()    
        done = False
        rewardSum = 0

        while not(done):              
            action = greedPolicy(qValues, state, explorationRate)
            next_state, reward, done = preformAction(action)
            rewardSum += reward   
            target = reward + 0.9 * np.max(qValues[next_state])
            error = target - qValues[state][action]
            qValues[state][action] += learnRate * error
            state = next_state

            if shouldDrawToGrid:
                drawToGrid(qValues, action=actions[action], shouldColor=True)

            
        episodeRewards.append(rewardSum)
    return episodeRewards, qValues

qRewards, qValues = qLearning(grid, gamma=0.9, learnRate=1, shouldDrawToGrid=True)
drawToGrid(qValues, shouldColor=True)
qRewards, _ = zip(*[qLearning(grid, shouldDrawToGrid=False, explorationRate=0.1, learnRate=1) for _ in range(10)])
np.mean(qRewards)

avgRewards = np.mean(qRewards, axis=0)
meanRewards = [np.mean(avgRewards)] * len(avgRewards)

fig, ax = plt.subplots()
ax.set_xlabel('Episodes')
ax.set_ylabel('Reward')
ax.plot(avgRewards)
print('Mean Reward: {}'.format(meanRewards[0]))

def main(qValues):

    state = reset()
    done = False

    while not(done):    
        action = greedPolicy(qValues, state, 0.0)
        next_state, reward, done = preformAction(action)      
        state = next_state  
        drawToGrid(qValues=qValues, action=actions[action], shouldColor=True)

main(qValues)

plt.pause(100)

