# nsockg

A simple script to build a large knowledge graph with biomedical relationships. Rebuild any time with:

```shell
python main.py
```

How it works:

This script is based on several tools:

1. [`bioversions`](https://github.com/cthoyt/bioversions) for automated discovery of new versions of data
2. [`pystow`](https://github.com/cthoyt/pystow) for automated downloading and caching of data
3. Years of experince in the chemistry lab, making me Not Scared of Chemistry
4. Years of experience with data processing and biological databases
5. [`zenodo-client`](https://github.com/zenodo-client) for automated uploading of the results to Zenodo, including automated versioning
