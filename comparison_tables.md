# Reference for xarray

## Comparison table
|function|xarray on python|xarray on MATLAB<br>(This package)|
| ---- | ---- | ---- |
| Create|DataArray(data) |DataArray(data)  |
|Lookup| (see below)|(see below)|
|Sum| x.sum('time') | x.sum('time') *|



## Lookup in DataArray
|Dimension lookup|Index lookup|python|MATLAB|
|----|----|----|----|
|Positional|By integer|da[:,0]|da(:,0)|
|Positional|By label|da[:,'IA']||
|By name|By integer|da.isel(space=0)<br>da[dict(space=0)]|da(space=1) *|
|By name|By label|da.sel(space='IA')<br>da.loc[dict(space='IA')]|da.sel(space='IA')<br>da.sel('space', 'IA')|
