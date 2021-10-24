# nsockg

A simple script to build a large knowledge graph with biomedical relationships.
Rebuild any time with `tox` using the following commands in your shell:

```shell
$ git clone https://github.com/cthoyt/nsockg.git
$ cd nsockg
$ pip install tox
$ tox
```

This script follows the [Bio2BEL philosophy](https://doi.org/10.1101/631812) of
completely automating data acquisition, storage, processing, and output of
biological data sources. It is based on several tools:

1. [`bioversions`](https://github.com/cthoyt/bioversions) for automated
   discovery of new versions of data
2. [`pystow`](https://github.com/cthoyt/pystow) for automated downloading and
   caching of data
3. [`zenodo-client`](https://github.com/zenodo-client) for automated uploading
   of the results to Zenodo, including automated versioning

## Download

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4597866.svg)](https://doi.org/10.5281/zenodo.4597866)
