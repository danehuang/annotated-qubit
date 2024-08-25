# The Annotated Qubit 

Web version: [(https://danehuang.github.io/annotated-qubit/)](https://danehuang.github.io/annotated-qubit/)

We provide an introduction to quantum computing (QC) and quantum information (QI).
The series of notes are partially based off of an advanced undergraduate/graduate
course taught at San Francisco State University (SFSU). The 
[series of notebooks](https://danehuang.github.io/annotated-qubit/) has
several parts:
1. Foundations: introduction to single qubit systems.
2. Foundations II: introduction to multi-qubit systems.
3. QI: basic quantum information.
4. QC: basic quantum computing algorithms.
5. Shor's: Shor's algorithm, i.e., advanced quantum computing.

Topic Dependency Graph
```
    Foundations
        |
  Foundations II
     /      \
    QC      QI
    |
  Shor's
```

## Installation

```
python -m venv annotated-qubit
source annotated-qubit/bin/activate
pip install -r requirements.txt
```


## Citation

If you found our work useful, please consider citing it:

```
@misc{huang2024annotatedqubit,
  author       = {Daniel Huang},
  title        = {The Annotated Qubit},
  year         = {2024},
  url          = {https://github.com/danehuang/annotated-qubit},
  note         = {GitHub repository},
}
```
