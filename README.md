# Quantum-Walk-on-Tree-Graphs
Code used for the publication: "Finding Paths in Tree Graphs With a Quantum Walk" (2018)

arXiv pre-print version: https://arxiv.org/abs/1710.05084

The code assmebled in this project were used in the paper mentioned above, which was the main focus of my thesis.  I've created this github project to share the code I created, in order to study these Scattering Quantum Random Walks.  For information about these quantum walks, I reccomended checking out the paper on arXiv.  All of the code in this project is designed for quantum random walks, specfificaly for the geometry of "N-Tree Graphs."

## Getting Started

All of the code is written for Python.  For an IDE to run the codes, I reccomend using Spyder through the Anaconda Distribution.

link: https://www.anaconda.com/download

The code requires Python 3.5 or higher

Once you have Spyder, or another Python IDE up and running, download all of the python files and put them together in a location somwhere on your computer. For example, store them in a folder on your desktop labeled "NTree Quantum Walks."  It is important that all of the python files be stored in the same location, as they call upon each other frequently to import functions.

For example, many of the codes will call upon the file 'Ntreefunctions.py' for functions:

```
import Ntreefunctions as nt
```
This Github comes with 4 python files that needed to be imported by other codes:

```
Ntreefunctions.py 
```
```
Quantum_Walk_Simulation.py 
```
```
Regression_Fits.py 
```
```
U_Eigenvalues.py 
```

As a good first test to make sure everything runs properly: open or copy into a new file, the python file named "Ntree_Test_Code".  If the file runs properly, a messege should print saying that all of the functions imported correctly.

## Classical Simulation of Quantum Systems

When designing new quantum algorithms, often times it is useful to run simulations of the behavior of quantum systems on a clssical computer.  This is precisely what all of these codes do: simulate the results one could expect from running a Quantum Random Walk on N-Tree graphs.  The advantage of simulating these walks classically is the ability to store information about the state of the system at any given moment.  Many of the codes in this project do exactly that: highlight the unique features of these quantum systems, with exact values for state amplitudes, probabilities, etc.

By studying the "under the hood" properties of these quantum systems, we can better determine whether they have the potential for speedups over classical algorithms.  Specifically, classical algorithms allow us to "freeze" these quantum systems at any point, and look at exact values for amplitudes / probabilities.  Such a task is impossible with real quantum systems, which is why studying them through classical codes is so insightful.

## Running The Codes
All of the codes provided in this project run "out of the box" and showcase certain properties of these N-Tree quantum systems.  Most of the codes produce a plot or print results to the terminal (or both).  For further explination on the results produced by individual codes, a short paragraph is provided at the beginning of each code.  For more information, I reccomended reading the arXiv paper listed above.  All of the plots from these coeds can be found in the paper, with better explinations.

## Coding Style Disclaimer

I am a physicist, not a professional software engineer.  If my codes fail to meet some coding etiquettes, I apologize!  I've spent quite a good deal of time making the codes as presentable and user-friendly as they are now, but I know they are still a little rough around the edges.


### Contact Me

**Daniel Koch** - dkochsjsu@gmail.com

If you have any questions / interests in the code, or quantum walks in general, feel free to reach out to me.
