# xarray_on_MATLAB
[![View xarray on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://jp.mathworks.com/matlabcentral/fileexchange/125995-xarray)

# Introduction
`xarray` is an N-dimensional labeled array, inspired by [xarray](https://github.com/pydata/xarray) implemented in Python.
The purpose of this `xarray` on MATLAB project is basically a port of the original (Python) xarray.

# Feature
Currently, `xarray` has a class of arrays called `DataArray`. Here is an example of the 2x3 2-dimensional `DataArray`
```
da = xarray.DataArray([10.1 10.4 12.3; 20.3 20.5 24.2] ...
, dims = {'x', 'time'}...
, coords = {'x', [2.5, 4.4], 'time', [0 .1 .2]})
```
## Accessing with Label
The most important feature of `DataArray` is the ability to select array elements by axis name. From the data above, we slice out the array for 
```
da.sel(x=4.4)
```
This selection method is revolutionary in two ways.
1. You do not have to remember the dimension number of the axis (in this case, `x` is the first dimension of the data array).
2. No need to calculate the index on the axis that is desired value (in this case, the second element of the `x` axis).

The second point could be realized by writing `matlab_array(:,x_vec==4.4)` and so on. However, not having to remember the dimension of the axis is useful when dealing with data that has many dimensions.
## Accessing as MATLAB style
The element of DataArray can be accessed also by conventional MATLAB array style. This is an example which has same result as written above.
```
data(2,:)
```
# Detailed Usage
See comparison_table.md for a table comparing how `xarray` on MATLAB is written with the original `xarray` notation.
Also, see readme.mlx for tutorial.

# Installation
Just copy +xarray folder to your folder contained the code.