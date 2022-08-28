import sys
import matplotlib.pyplot as plt
import Lablibrary

sys.stdin = open("input.txt", "r")
sys.stdout = open("output.txt", "w")
sys.stderr = open("error.txt", "w")

# question 2
seed1, seed2, seed3, N = list(map(float,input().strip().split()))[:4]
vol = Lablibrary.volapprox(seed1, seed2, seed3, N)
print(vol)

# question 3
seed1, seed2, N = list(map(float,input().strip().split()))[:3]

print('RMS                   Discplacement')

x_value, y_value, rms, disc = Lablibrary.RandomWalk(seed1, seed2, N)
Lablibrary.plot(x_value, y_value,'RandomWalk_300')
print('{}, {}'.format(rms, disc))

seed1, seed2, N = list(map(float,input().strip().split()))[:3]

x_value, y_value, rms, disc= Lablibrary.RandomWalk(seed1, seed2, N)
Lablibrary.plot(x_value, y_value, 'RandomWalk_600')
print('{}, {}'.format(rms, disc))

seed1, seed2, N = list(map(float,input().strip().split()))[:3]

x_value, y_value, rms, disc = Lablibrary.RandomWalk(seed1, seed2, N)
Lablibrary.plot(x_value, y_value,'RandomWalk_900')
print('{}, {}'.format(rms, disc))
