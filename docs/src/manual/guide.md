# Quick Guide
If what goes inside Li-ion batteries interests you and you want to simulate various possible parameters for a given input to a cell, this package is for you!

## Installation
`LiionBatteryModels.jl` can be installed in 2 ways as below:

### From Julia package manager

The simplest way to install is to install the package from the julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> add CaseStudy
```

### From Github 
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> add https://github.com/just-ary27/LiionBatteryModels.jl.git
```

## Using the models
Presently there are 2 models implemented in the library. These are as below:

---

### SenthilModel
This is a reduced order model for a lithium ion cell with uniform reaction rate approximation. To start using the model, import it from the CaseStudy core. We abbreviate it as `sm` for convinience.

```julia
import CaseStudy.SenthilModel as sm
```

### AshwiniModel
This is a closed form reduced order electrochemical model for lithium ion cells. To start using it import it from CaseStudy core like previously. We abbreviate it as `am` for convinience.

```julia
import CaseStudy.SenthilModel as am
```

---

*For a detailed guide on how to use each model, see [examples.](../manual/examples.md)*