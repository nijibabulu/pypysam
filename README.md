A python reimplementation of htslib functions, making it compatible with `pypy`. Parsing and writing does not perform as well as its htslib counterparts. Limited documentation available, please see the readthedocs pages for pysam:

http://pysam.readthedocs.io/en/latest/

To install, use 

```pypy setup.py install```

If you would like to install a local copy in `/home` rather than `/usr` space, use
`pypy setup.py install --prefix=$HOME/.local` and update your `$PYTHONPATH` variable
to include `$HOME/.local/lib/pythonX.X/site-packages`.

See `pypy setup.py --help` for further instructions.

