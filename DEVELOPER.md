## Notes for Dev's

### To upload to PyPI

Currently running in the main `nuctools` dir:
```
python -m pip install --upgrade build  
python -m build
twine check dist/nuctools-<version>
twine upload dist/nuctools-<version>
```

and then using the API Token method for authentication. 

