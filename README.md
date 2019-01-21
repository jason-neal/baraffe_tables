# Baraffe tables.
[![Build Status](https://travis-ci.org/jason-neal/baraffe_tables.svg?branch=master)](https://travis-ci.org/jason-neal/baraffe_tables)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/cdf265c880db4d3289adf77c41d5f007)](https://www.codacy.com/app/jason-neal/baraffe_tables?utm_source=github.com&utm_medium=referral&utm_content=jason-neal/baraffe_tables&utm_campaign=badger)[![Updates](https://pyup.io/repos/github/jason-neal/baraffe_tables/shield.svg)](https://pyup.io/repos/github/jason-neal/baraffe_tables/)[![Python 3](https://pyup.io/repos/github/jason-neal/baraffe_tables/python-3-shield.svg)](https://pyup.io/repos/github/jason-neal/baraffe_tables/)

`baraffe_tables` is a small package to access the Baraffe evolutionary tables and interpolating between ages, and rows.
It also calculates the flux ratio of suspected companions given their mass.

How to use
----------
To install:
```python
python setup.py install
```
from the cloned repo.

To access a parameter from the Baraffe tables
```python
from baraffe_tables.table_search import baraffe_table_search
row = baraffe_table_search(column="M/Ms", value=0.08, age=4.0, model=2003, age_interp=True)
print(row)
```
or you can use the `query_baraffe.py` CLI script. Append `-h`for help.


Host-Companion flux ratio
-------------------------
There are also two scripts to use the [Baraffe et al. 2003](http://adsabs.harvard.edu/abs/2003A%26A...402..701B)/[2015](http://adsabs.harvard.edu/abs/2015A%26A...577A..42B) to estimate Brown Dwarf (BD) parameters knowing some information...
These are not fixed to only BD compaions but can be planetary or stellar companions also, withing the Baraffe models mass range.

- `mass_to_flux_ratio.py`

    Takes a companion mass (in Mjup) and uses this to look up the absolute magnitude of the brown dwarf. This compares it to the magnitude of the star in the literature to calculate the flux ratio.

- `flux_ratio_to_mass.py`

    This attempts to go back the opposite direction. Given  the flux ratio what is the mass of the companion.
