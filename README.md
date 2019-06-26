# motifHeat 

Wrapper to format heamap from defined [GNPS-MZmine](https://gnps.ucsd.edu/) inputs

### Table of contents

* [Installation](#installation)
* [Run motifHeat](#main_citation)
* [License](#license)

## Installation 

Install motifHeat with:

 `devtools::install_github("computational-chemical-biology/motifHeat")`

## Run motifHeat 

To run `motifHeat` simply do:

```
h <- format_heatmap(tab, meta, selectField='StrainName', selectValue='Burkholderia dolosa AU0645  Genomovar type VI', factorColList=factorColList)
```

You should be able to generate the following image

<img src="img/example_heatmap.png"/>

## License
This repository is available under the following license https://github.com/
