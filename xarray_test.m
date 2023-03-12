sx = 2; sy = 10; sz = 3; st = 5;
temperature = reshape(1:(sx*sy*sz*st), [sx, sy, sz, st]);
lon = reshape(2*(1:(sx*sy*sz)), [sx, sy, sz]);
lat = reshape(3*(1:(sx*sy*sz)), [sx, sy, sz]);
time = 4*(1:st);
x = 5*(1:sx);
y = 6*(1:sy);
z = 7*(1:sz);
reference_time = 2014;
%% coordinates with dimensions
da = xarray.DataArray(temperature,...
    dims = {'x', 'y','z','time'},...
    coords = {'lon', {{'x','y','z'},lon}, 'lat', {{'x','y','z'}, lat}, 'time', time, 'reference_time', reference_time}, ...
    attrs = {'description', 'Ambient temerature.','units','degC'});

testFileName = 'test/da1';
refFileName = [testFileName, '_ref'];
if(exist(testFileName, 'file'))
    delete(testFileName);
end
diary(testFileName)
disp(da);
diary off
flgIdentical = diffFile(testFileName, refFileName);
if(~flgIdentical)
    visdiff(refFileName, testFileName);
end
assert(flgIdentical);
%% not specify dimension name
da = xarray.DataArray(temperature);

testFileName = 'test/da2';
refFileName = [testFileName, '_ref'];
if(exist(testFileName, 'file'))
    delete(testFileName);
end
diary(testFileName)
disp(da);
diary off
flgIdentical = diffFile(testFileName, refFileName);
if(~flgIdentical)
    visdiff(refFileName, testFileName);
end
assert(flgIdentical);

%% lack of dim's names
catchError = false;
try
    da = xarray.DataArray(temperature,...
        dims = {'x', 'y','z'},...%need 4 dimensions
        coords = {'lon', {{'x','y','z'},lon}, 'lat', {{'x','y','z'}, lat}, 'time', time, 'reference_time', reference_time}, ...
        attrs = {'description', 'Ambient temerature.','units','degC'});
catch MExc
    catchError = true;
    assert(MExc.identifier == "xarray:invaliddimensions");
end
assert(catchError);

%% too many dim's names -> Do not treat as error
da = xarray.DataArray(temperature,...
    dims = {'x', 'y','z','time', 'w'},...%need 4 dimensions
    coords = {'lon', {{'x','y','z'},lon}, 'lat', {{'x','y','z'}, lat}, 'time', time, 'reference_time', reference_time}, ...
    attrs = {'description', 'Ambient temerature.','units','degC'});


%% wrong dimension name in coordinate list
catchError = false;
wrong_lon = reshape(2*(1:(sx*(sy-1)*sz)), [sx, (sy-1), sz]);

try
da = xarray.DataArray(temperature,...
    dims = {'x', 'y','z','time'},...
    coords = {'lon', {{'x','y','z'},wrong_lon}, ...%wrong dimension
    'lat', {{'x','y','z'}, lat}, 'time', time, 'reference_time', reference_time}, ...
    attrs = {'description', 'Ambient temerature.','units','degC'});
catch MExc
    catchError = true;
    assert(MExc.identifier == "xarray:invalidcoorddims");
end
assert(catchError);

%% wrong dimension in coordinate list
catchError = false;
try
da = xarray.DataArray(temperature,...
    dims = {'x', 'y','z','time'},...
    coords = {'lon', {{'x','y','w'},lon}, ...%wrong dimension name
    'lat', {{'x','y','z'}, lat}, 'time', time, 'reference_time', reference_time}, ...
    attrs = {'description', 'Ambient temerature.','units','degC'});
catch MExc
    catchError = true;
    assert(MExc.identifier == "xarray:invaliddimensionname");
end
assert(catchError);

%% 2 dimensional empty
da = xarray.DataArray(double.empty(0,0));

testFileName = 'test/empty2';
refFileName = [testFileName, '_ref'];
if(exist(testFileName, 'file'))
    delete(testFileName);
end
diary(testFileName)
disp(da);
diary off
flgIdentical = diffFile(testFileName, refFileName);
if(~flgIdentical)
    visdiff(refFileName, testFileName);
end
assert(flgIdentical);
assert(isscalar(da)==false);

da = xarray.DataArray.empty(0,0);

testFileName = 'test/empty2';
refFileName = [testFileName, '_ref'];
if(exist(testFileName, 'file'))
    delete(testFileName);
end
diary(testFileName)
disp(da);
diary off
flgIdentical = diffFile(testFileName, refFileName);
if(~flgIdentical)
    visdiff(refFileName, testFileName);
end
assert(flgIdentical);
assert(isscalar(da)==false);

da = xarray.DataArray.empty;

testFileName = 'test/empty2';
refFileName = [testFileName, '_ref'];
if(exist(testFileName, 'file'))
    delete(testFileName);
end
diary(testFileName)
disp(da);
diary off
flgIdentical = diffFile(testFileName, refFileName);
if(~flgIdentical)
    visdiff(refFileName, testFileName);
end
assert(flgIdentical);
assert(isscalar(da)==false);


%% 1 dimensional empty
da = xarray.DataArray(double.empty(0), dims={'x'});

testFileName = 'test/empty1';
refFileName = [testFileName, '_ref'];
if(exist(testFileName, 'file'))
    delete(testFileName);
end
diary(testFileName)
disp(da);
diary off
flgIdentical = diffFile(testFileName, refFileName);
if(~flgIdentical)
    visdiff(refFileName, testFileName);
end
assert(flgIdentical);
assert(isscalar(da)==false);

%% 3 dimensional empty
da = xarray.DataArray.empty(0,1,0);

testFileName = 'test/empty3(0_1_0)';
refFileName = [testFileName, '_ref'];
if(exist(testFileName, 'file'))
    delete(testFileName);
end
diary(testFileName)
disp(da);
diary off
flgIdentical = diffFile(testFileName, refFileName);
if(~flgIdentical)
    visdiff(refFileName, testFileName);
end
assert(flgIdentical);
assert(isscalar(da)==false);

da = xarray.DataArray.empty(1,0,0);

testFileName = 'test/empty3(1_0_0)';
refFileName = [testFileName, '_ref'];
if(exist(testFileName, 'file'))
    delete(testFileName);
end
diary(testFileName)
disp(da);
diary off
flgIdentical = diffFile(testFileName, refFileName);
if(~flgIdentical)
    visdiff(refFileName, testFileName);
end
assert(flgIdentical);
assert(isscalar(da)==false);

%% one dimensional
da = xarray.DataArray(2:2:6, dims={'x'},coords={'x', 0.5:0.5:1.5});

testFileName = 'test/onedim';
refFileName = [testFileName, '_ref'];
if(exist(testFileName, 'file'))
    delete(testFileName);
end
diary(testFileName)
disp(da);
disp(da(2:3));
da_trim = da([1,3]);
disp(da_trim);
da(2) = [];
disp(da);
diary off
flgIdentical = diffFile(testFileName, refFileName);
if(~flgIdentical)
    visdiff(refFileName, testFileName);
end
assert(flgIdentical);
assert(isscalar(da)==false);

%% one dimensional one element
da = xarray.DataArray(2, dims={'x'},coords={'x', 0.5});

testFileName = 'test/onedimoneel';
refFileName = [testFileName, '_ref'];
if(exist(testFileName, 'file'))
    delete(testFileName);
end
diary(testFileName)
disp(da);
diary off
flgIdentical = diffFile(testFileName, refFileName);
if(~flgIdentical)
    visdiff(refFileName, testFileName);
end
assert(flgIdentical);
assert(isscalar(da)==true);

%% two dimensional one element
da = xarray.DataArray(2, dims={'x','y'},coords={'x', 0.5});

testFileName = 'test/twodimoneel';
refFileName = [testFileName, '_ref'];
if(exist(testFileName, 'file'))
    delete(testFileName);
end
diary(testFileName)
disp(da);
diary off
flgIdentical = diffFile(testFileName, refFileName);
if(~flgIdentical)
    visdiff(refFileName, testFileName);
end
assert(flgIdentical);
assert(isscalar(da)==true);

%% three dimensional one element
da = xarray.DataArray(2*ones(1,1,1), dims={'x','y','z'},coords={'x', 0.5});

testFileName = 'test/3Doneel';
refFileName = [testFileName, '_ref'];
if(exist(testFileName, 'file'))
    delete(testFileName);
end
diary(testFileName)
disp(da);
diary off
flgIdentical = diffFile(testFileName, refFileName);
if(~flgIdentical)
    visdiff(refFileName, testFileName);
end
assert(flgIdentical);
assert(isscalar(da)==true);


%% matrix multiple scalar
da = xarray.DataArray([1 2;3 4]);
da2 = xarray.DataArray(3);
mult = da*3;
assert( all(mult.data(:).'==[1 3 2 4]*3));
mult = 3*da;
assert( all(mult.data(:).'==[1 3 2 4]*3));
mult = da*da2;
assert( all(mult.data(:).'==[1 3 2 4]*3));
mult = da2*da;
assert( all(mult.data(:).'==[1 3 2 4]*3));

%% matrix multiple matrix
da = xarray.DataArray([1 2 3;3 4 5]);
da2 = xarray.DataArray([1 2;3 4]);
answer = [1 2;3 4]*[1 2 3;3 4 5];
mult = [1 2;3 4]*da;
assert(all(mult.data==answer,'all'))
mult = da2*da;
assert(all(mult.data==answer,'all'))

da = xarray.DataArray([1 2 ;3 3; 4 5]);
da2 = xarray.DataArray([1 2;3 4]);
answer = [1 2 ;3 3; 4 5]*[1 2;3 4];
mult = da*[1 2;3 4];
assert(all(mult.data==answer,'all'))
mult = da*da2;
assert(all(mult.data==answer,'all'))

%% multiple 2d matrix
da = xarray.DataArray([1 2 3;3 4 5]);
da2 = xarray.DataArray([1 2 3;3 4 5]);
answer = [1 2 3;3 4 5].*[1 2 3;3 4 5];
mult = da.*da2;
assert(all(mult.data==answer,'all'))

%% slice
da = xarray.DataArray(temperature,...
    dims = {'x', 'y','z','time'},...
    coords = {'lon', {{'x','y','z'},lon}, 'lat', {{'x','y','z'}, lat}, 'time', time, 'reference_time', reference_time}, ...
    attrs = {'description', 'Ambient temerature.','units','degC'});

testFileName = 'test/slice';
refFileName = [testFileName, '_ref'];
if(exist(testFileName, 'file'))
    delete(testFileName);
end
diary(testFileName)
disp(da(:,[3,6:10], 1, 1:2:5));
diary off
flgIdentical = diffFile(testFileName, refFileName);
if(~flgIdentical)
    visdiff(refFileName, testFileName);
end
assert(flgIdentical);


%% delete elements
da = xarray.DataArray(temperature,...
    dims = {'x', 'y','z','time'},...
    coords = {'lon', {{'x','y','z'},lon}, 'lat', {{'x','y','z'}, lat}, 'time', time, 'reference_time', reference_time}, ...
    attrs = {'description', 'Ambient temerature.','units','degC'});
testFileName = 'test/delete';
refFileName = [testFileName, '_ref'];
if(exist(testFileName, 'file'))
    delete(testFileName);
end
diary(testFileName)
da(:,[2 4],:,:) = [];
disp(da);
diary off
flgIdentical = diffFile(testFileName, refFileName);
if(~flgIdentical)
    visdiff(refFileName, testFileName);
end
assert(flgIdentical);

%% append new value
da = xarray.DataArray([1 2;3 4], dims={'x','y'}, coords = {'x', [3.4, 3.5], 'y',[1.1, 1.4]});
da(3,3,2) = 1;
testFileName = 'test/append';
refFileName = [testFileName, '_ref'];
if(exist(testFileName, 'file'))
    delete(testFileName);
end
diary(testFileName)
disp(da);
diary off
flgIdentical = diffFile(testFileName, refFileName);
if(~flgIdentical)
    visdiff(refFileName, testFileName);
end
assert(flgIdentical);


function tf = diffFile(file1, file2)
    fid1 = fopen(file1,"r");
    assert(fid1 ~= -1,'file1 cannot open');
    fid2 = fopen(file2,"r");
    if(fid2 == -1)
        fclose(fid1);
        error('file2 cannot open');
    end
    try
        lineNo = 0;
        while(true)
            lineNo = lineNo +1;
            tline1 = fgets(fid1);
            tline2 = fgets(fid2);
            if(isnumeric(tline1) ~= isnumeric(tline2))
                tf = false;
                break;
            elseif(isnumeric(tline1))
                tf = true;
                break;
            end
            if(~strcmp(tline1, tline2))
                tf = false;
                break;
            end
        end
        if(~tf)
            fprintf('Different in line %d in file %s\n', lineNo, file1)
        end
    catch ME
        fclose(fid1);
        fclose(fid2);
        rethrow(ME);
    end
    fclose(fid1);
    fclose(fid2);
end