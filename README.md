# OEPandas Schrodinger Support

**Author:** Scott Arne Johnson ([scott.arne.johnson@gmail.com](mailto:scott.arne.johnson@gmail.com))

Schrodinger file format support for OEPandas

## Description

Coming soon!

## Development

This package uses ```pyproject.toml``` with ```setuptools``` to manage the package installation. All package
requirements should go in ```requirements.txt```.

You can install this package in "editable mode," which allows you to make it available in Python while you are
developing it. This facilitates both development and testing. To make editable mode compatible with most development
IDEs like PyCharm and Visual Studio Code, you can execute:

```bash
pip install --config-settings editable_mode=compat -e ".[dev]"
```

If you have the ```invoke``` package installed, you can automate deployment of your package (wheel and source
distributions) to the BMS pypi server with:

```bash
invoke upload
```