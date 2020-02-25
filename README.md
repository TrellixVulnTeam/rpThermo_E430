# rpThermodynamics

Galaxy tool that reads a single or collection of rpSBML files, parses the SMILES or InChI structures to calculate the formation energy of the molecules. Either calculates the Gibbs free energy using databases derived from [Equilibrator](https://gitlab.com/equilibrator/component-contribution) or [Alberty et al.](https://www.biorxiv.org/content/biorxiv/early/2018/12/12/492819.full.pdf). Since the retosynthesis algorithm generates chemical species without corss-references to databases, the tool falls back to the [Component Contribution](https://github.com/eladnoor/component-contribution) method to estimate the Gibbs free energy of a reaction.

## Information Flow

### Input

Required information:
* **Input rpSBML**: Either tar.xz collection of rpSBML or single rpSBML file

Advanced options:
* **Name of the heterologous pathway**: (default: rp_pathway) Groups ID of the heterologous pathway
* **IP address of the rpFBA REST service**: IP address of the REST service

### Output

* **rpThermo**: The output is a tar.xz archive containing a list of rpSBML files or a single SBML file

## Installing

To build the image using the Dockerfile, use the following image:

```
docker build -t brsynth/rpthermo-rest:dev .
```

To run the service as localhost use the following service:

```
docker run -p 8882:8888 brsynth/rpthermo-rest:dev
```

## Algorithm

The thermodyamics is calculated using the [Component Contribution](https://gitlab.com/equilibrator/component-contribution) method

## Prerequisites

* Docker - [Install](https://docs.docker.com/v17.09/engine/installation/)
* libSBML - [Anaconda library](https://anaconda.org/SBMLTeam/python-libsbml)
* Component Contribution - [Git to the project](https://gitlab.com/elad.noor/component-contribution)

## Contributing

TODO

## Versioning

Version 0.1

## Authors

* **Melchior du Lac** 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thomas Duigou
* Joan Hérisson

### How to cite rpThermodynamics?

TODO
