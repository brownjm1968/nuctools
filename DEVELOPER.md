## Notes for Dev's

### To contribute: 

1. Create a fork
2. Make changes to your forked repo
3. Create pull request to merge into the main repo here

### Documentation

All documentation should be in the numpydoc style for functions and classes 
or written into the RST format files in `docs/`

### Testing

Any new functions should come with an accompanying test in the `tests/` folder
with file name matching `test_*.py`. 

### To manually upload to PyPI

Currently running in the main `nuctools` dir:
```
python -m pip install --upgrade build  
python -m build
twine check dist/nuctools-<version>
twine upload dist/nuctools-<version>
```

and then using the API Token method for authentication. 

