# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 00:49:12 2024

@author: comel
"""

import matplotlib.pyplot as plt

eigenvalues = []

with open('eigenvalues.txt', 'r') as file:
    for line in file:
        eigenvalues.append(float(line.strip()))

index = list(range(len(eigenvalues)))

plt.scatter(index, eigenvalues)
plt.xlabel("Index")
plt.ylabel("Eigenvalues")
plt.title("Plot of eigenvalues")

plt.show()
