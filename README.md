# Quantum-Walk-on-Tree-Graphs
Code used for the publication: "Finding Paths in Tree Graphs With a Quantum Walk" (2018)

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

As a good first test to make sure everything runs properly: open or copy into a new file, the python file named "Ntree Test Code".  If the file runs properly, a messege should print saying that all of the functions imported correctly.

## Classical Simulation of Quantum Systems

When designing new quantum algorithms, often times it is useful to run simulations of the behavior of quantum systems on a clssical computer.  This is precisely what all of these codes do: simulate the results one could expect from running a Quantum Random Walk on N-Tree graphs.  The advantage of simulating these walks classically is the ability to store information about the state of the system at any given moment.  Many of the codes in this project do exactly that: highlight the unique features of these quantum systems, with exact values for state amplitudes, probabilities, etc.

By studying the "under the hood" properties of these quantum systems, we can better determine whether they have the potential for speedups over classical algorithms.

--------------------------


## Reccomended Files to Run

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc
