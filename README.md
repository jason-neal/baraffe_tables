# Baraffe tables.
[![Build Status](https://travis-ci.org/jason-neal/baraffe_tables.svg?branch=master)](https://travis-ci.org/jason-neal/baraffe_tables)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/cdf265c880db4d3289adf77c41d5f007)](https://www.codacy.com/app/jason-neal/baraffe_tables?utm_source=github.com&utm_medium=referral&utm_content=jason-neal/baraffe_tables&utm_campaign=badger)[![Updates](https://pyup.io/repos/github/jason-neal/baraffe_tables/shield.svg)](https://pyup.io/repos/github/jason-neal/baraffe_tables/)[![Python 3](https://pyup.io/repos/github/jason-neal/baraffe_tables/python-3-shield.svg)](https://pyup.io/repos/github/jason-neal/baraffe_tables/)

`baraffe_tables` is a small package to access the Baraffe (2003, 2015) evolutionary tables and enables interpolating between rows of the tables and between stellar ages (different tables).
It was primarily created to calcualte the flux ratio of stellar companions given their mass.

How to use
----------
To clone and install:
```bash
git clone https://github.com/jason-neal/baraffe_tables
cd baraffe_tables
python setup.py install
```

`baraffe_table_search` can be used to access a particular row of the Baraffe tables given a value of one parameter.
For example, the parameters of a star with a mass of M/Ms=0.08 and an age=4.0 can be obtained using:
```python
from baraffe_tables.table_search import baraffe_table_search
row = baraffe_table_search(column="M/Ms", value=0.08, age=4.0, model=2003, age_interp=True)
print(row)
```
The `query_baraffe.py` CLI script can also be used to achieve the same. Append `-h`for help.
```bash
query_baraffe.py -h
```


Host-Companion flux ratio
-------------------------
There are also two scripts to use the [Baraffe et al. 2003](http://adsabs.harvard.edu/abs/2003A%26A...402..701B)/[2015](http://adsabs.harvard.edu/abs/2015A%26A...577A..42B) to estimate stellar companion parameters (primarily Brown Dwarfs (BD)) knowing their mass.
These can be used for any companion within the Baraffe models mass range.

- `mass_to_flux_ratio.py`

    Takes a companion mass (in Mjup) and uses this to look up the absolute magnitude of the companion. This compares the magnitude of companion to the magnitude of the star in the literature (SIMBAD) to calculate the flux ratio. This also accounts for the conversion between apparent and absolute magnitudes.

- `flux_ratio_to_mass.py`

    This script works in reverse to  `mass_to_flux_ratio.py`. Given a companion/host flux ratio what is the mass of the companion.

For example:
HD30501 has a BD companion with a mass of 89 Mjup and an age of 5Gyr. The companion/host flux ratio an be estimated using the 2003 models with:
```bash
mass_to_flux_ratio.py HD30501 89 5 -m 03
```

Contributing
-------------
Any issues or suggestions?
Everyone is welcome to make an [issue](https://github.com/jason-neal/baraffe_tables/issues) or a [pull request](https://github.com/jason-neal/baraffe_tables/pulls). 
