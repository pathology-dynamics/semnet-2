
# SemNet

This repository provides source code, tests, and examples for SemNet 2.0. 

SemNet (now SemNet 2.0) is a group of modules for working with semantic networks, specifically those built using the Semantic MEDLINE Database (SemMedDB) repository of semantic predications. These predications (subject-predicate-object triples) are modeled as nodes (subject and object) and edges (predicate) to create an expansive biomedical concept graph. Algorithms (detailed in the SemNet 2.0 paper) were then implemented to rank nodes in terms of their connections to other nodes, allowing for intuitive and novel literature based discovery. 

Why SemNet 2.0?
* Improvements on algorithm efficiency.
* Re-engineered graph data structure.
* Replacement of ULARA and improved rank aggregation.
* Potential 1,000,000,000 fold speed improvement compared to the original SemNet (see SemNet 2.0 paper for more precisely quantified results).

## Table of Contents

* [Technologies](#technologies)
* [Installation](#installation)
* [Usage](#usage)
* [Project Status](#project-status)
* [Contributing](#contributing)
* [Credits](#credits)
* [License](#license)

## Technologies

SemNet 2.0 uses the following technologies:

* Python == 3.9.6
  * numpy == 1.20.3

The examples and testing scripts use the following technologies:

* Python == 3.9.6
  * numpy == 1.20.3
  * Pandas == 1.3.3
  * lxml == 4.6.3
  * requests == 2.26.0

## Installation

1. To install SemNet 2.0, clone the repository to your machine.
2. Change to the semnet directory.
3. Use the command `pip install .` (include the period).
4. You should now be able to `import semnet` and use the package.

## Usage

To see an example notebook using SemNet 2.0, navigate to `/examples/biomedical_kg_example.ipynb`.

Unlike the original SemNet, SemNet 2.0 is designed to work with predications (derived from SemMedDB) specifically structured as a list of 'records' formatted Python dictionaries (see: https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_dict.html). `/examples/sample_data.csv` shows the first ten rows of a larger .csv file used in-house that has all of the necessary columns for the HetGraph object to create the knowledge graph (in this particular .csv each row corresponds to a predication). Converting the .csv file into the list of dictionaries can also be observed in `/examples/biomedical_kg_example.ipynb`. Finally, an inverse relationship dictionary is required for SemNet 2.0. See `rel2inv` in `biomedical_kg_example.ipynb`.

Useful Links
* SemNet 2.0 paper: tbd
* Original SemNet paper: <https://www.frontiersin.org/articles/10.3389/fbioe.2019.00156/full>
* HeteSim paper: <https://arxiv.org/abs/1309.7393>
* SemRep and SemMedDB: <https://lhncbc.nlm.nih.gov/ii/tools/SemRep_SemMedDB_SKR.html>

## Project Status

This project is currently released for public use and is being actively maintained.

## Contributing

If you feel that this repo can be improved in any way and would like to contribute, contact Anna Kirkpatrick <akirkpatrick3@gatech.edu>, Stephen Allegri <sallegri3@gatech.edu>, and Cassie Mitchell <cassie.mitchell@bme.gatech.edu>.

## Credits

To cite this project and/or code, please use `citation.bib`, which is provided in the repository. 

Alternatively, use GitHub's built in citation feature / CITATION.cff.

## License

This project uses the GNU AGPLv3 license. In addition to the terms of this license, any work using this code must cite this repository using the instructions above.
