classdef DataArray < matlab.mixin.indexing.RedefinesParen
    properties (Access=public)
        cont_array % 多次元配列．
        sz % 整数ベクトル．cont_arrayの各次元の個数としてユーザー側が認識しているもの
        % 基本的にはsize(cont_array)．ただし，Matlabでは1次元のcont_arrayが存在しないため，1次元を指定された場合のみ，size(cont_array)と異なる．
        ndims % 整数. cont_arrayの次元の数としてユーザー側が認識しているもの．= numel(sz)
        ncoords % 整数. 定義されている軸の数．原理上いくらでも．
        dims_name % string配列．dimensionの名前. numel(dims_name)=ndims
        cid % containers.Map 軸名→coords ID
        coords % セル配列．軸の値を与える配列
        coords_name % string配列．
        coords_isSpecified % logical配列．coordsに具体的な値が指定されているかどうか．
        coords_dims % 2次元logical配列．size(coords_dims) = [ncoords, ndims]
        attrs
    end
    
    methods
        function obj = DataArray(data, NameValueArgs)
            arguments
                data = []
                NameValueArgs.dims = {}
                NameValueArgs.coords = {}
                NameValueArgs.attrs = {}
            end
            argDims = NameValueArgs.dims;
            argCoords = NameValueArgs.coords;
            argAttrs = NameValueArgs.attrs;

            obj.attrs = argAttrs;
            
            % Setting on contained array
            % Setting on dimensions
            temp_sz = size(data);
            temp_ndims = numel(temp_sz);
            obj.dims_name = string.empty;
            if(numel(argDims) == 0) %dimension is determined by data.
                for idxDims = 1:temp_ndims
                    obj.dims_name(idxDims) = sprintf("dims_%d",idxDims);
                end
                obj.cont_array = data;
                obj.ndims = temp_ndims;
            else
                obj.dims_name = string(argDims);
                for idxDims = 1:numel(obj.dims_name)
                    if(strlength(obj.dims_name(idxDims)) == 0)
                        obj.dims_name(idxDims) = sprintf("dims_%d",idxDims);
                    end
                end
                if(temp_ndims == 2 && numel(obj.dims_name) == 1)
                    %Matlabが1次元配列を使えないことによる特殊ケースとして取り扱う
                    obj.cont_array = data(:); %縦ベクトルにしておく
                    obj.ndims = 1;
                else 
                    assert(temp_ndims <= numel(obj.dims_name)...
                        ,'xarray:invaliddimensions'...
                        ,"different number of dimensions on data and specified in dims %d vs %d", temp_ndims, numel(obj.dims_name));
                    obj.cont_array = data;
                    obj.ndims = numel(obj.dims_name);
                end
            end
            obj.sz = obj.csize();

            % Setting on coordinates
            % set dimension as the first coordinates
            obj.cid = containers.Map;
            obj.coords = {};
            obj.coords_name = string.empty;
            obj.coords_isSpecified = logical.empty;
            obj.coords_dims = logical.empty;
            for idxDim = 1:obj.ndims
                obj.cid(obj.dims_name(idxDim)) = idxDim;
                obj.coords_name(idxDim) = obj.dims_name(idxDim);
                %Give an integer sequence as coodinates to the axes that are the Dimension at this moment.
                obj.coords_isSpecified(idxDim) = false;
                obj.coords{idxDim} = transpose(1:obj.sz(idxDim));
                obj.coords_dims(idxDim, idxDim) = true;
            end
            obj.ncoords = obj.ndims;
            for idxArgCoord = 1:2:numel(argCoords)
                r_coord_name = argCoords{idxArgCoord};
                r_coord_info = argCoords{idxArgCoord+1};
                if(any(strcmp(obj.dims_name, r_coord_name))) %dimensionでもある軸の場合
                    tmp = find(strcmp(obj.dims_name, r_coord_name));
                    idxDim = tmp(1);
                    obj.coords_isSpecified(idxDim) = true;
                    obj.coords{idxDim} = r_coord_info(:);%縦ベクトルにする
                else % dimensionではないが，軸を与えられている場合．
                    obj.ncoords = obj.ncoords+1;
                    obj.cid(r_coord_name) = obj.ncoords;
                    obj.coords_name(obj.ncoords) = r_coord_name;
                    obj.coords_isSpecified(obj.ncoords) = true;
                    if(iscell(r_coord_info)) %軸としてdimensionを参照している場合．
                        r_dim_names = r_coord_info{1};
                        for idxRdim = 1:numel(r_dim_names)
                            assert(isKey(obj.cid, r_dim_names{idxRdim}), ...
                                'xarray:invaliddimensionname',...
                                '"%s", specified in coord:"%s", is not included in the dimension lists.', r_dim_names{idxRdim}, r_coord_name);
                            obj.coords_dims(obj.ncoords,obj.cid(r_dim_names{idxRdim})) = true;
                        end
                        if(numel(r_dim_names) == 1)
                            %指定した次元が1の場合は縦ベクトルにしておく．
                            obj.coords{obj.ncoords} = reshape(r_coord_info{2},[],1);
                        else
                            obj.coords{obj.ncoords} = r_coord_info{2};
                        end
                    else
                        % dimsが指定されてない場合，coordsがどのような形で与えられていようと問題ない．
                        obj.coords{obj.ncoords} = r_coord_info;
                        obj.coords_dims(obj.ncoords,:) = false; % すべてfalseの行を増やす
                    end
                end
            end
            obj.CheckConsistency();
        end
        
        function obj = expand_dims(obj, NameValueArgs)
            arguments
                obj
                NameValueArgs.dims = {};
                NameValueArgs.axis = obj.ndims+1;
            end
            argDims = NameValueArgs.dims;

            new_coords_name = string.empty;
            if(numel(argDims) == 0)
                nadddims = 1;
                newcid = NameValueArgs.axis;
                new_coords = {1};
                new_coords_name(1) = '';
                new_coords_isSpecified = false;
            else
                nadddims = numel(argDims);
                newcid = NameValueArgs.axis:NameValueArgs.axis+nadddims-1;
                new_coords = cell(0);
                new_coords_isSpecified = logical.empty;
                for idxAddDim = 1:nadddims
                    if(iscell(argDims{idxAddDim}))
                        new_coords_name(idxAddDim) = argDims{idxAddDim}{1};
                        new_coords(idxAddDim) = argDims{idxAddDim}(2);
                        new_coords_isSpecified(idxAddDim) = true;
                    else
                        new_coords_name(idxAddDim) = string(argDims);
                        new_coords(idxAddDim) = {1};
                        new_coords_isSpecified(idxAddDim) = false;
                    end
                end
            end
            for idxAddDim = 1:(nadddims-1)
                obj = obj.expand_single_dim(new_coords_name(idxAddDim), newcid(idxAddDim) ...
                    , new_coords{idxAddDim}, new_coords_isSpecified(idxAddDim));
                % new size setting and expand contained array
                reshapesize = [obj.sz(1:newcid-1), 1, obj.sz(newcid:end)];
                obj.cont_array = reshape(obj.cont_array, reshapesize);
                newsize = [obj.sz(1:newcid-1), numel(coord), obj.sz(newcid:end)];
                obj.cont_array = xarray.DataArray.extendArray(obj.cont_array, num2cell(newsize));
            end

        end
        
        function obj = expand_single_dim(obj, dimName, newcid, coord, isSpecified)
            if(strlength(dimName)==0) %dimension is determined by data.
                dimName = sprintf("dims_%d",newcid);
                cnter = 1;
                while(isKey(obj.cid, dimName))
                    dimName = sprintf("dims_%d",newcid+cnter-1);
                    cnter = cnter+1;
                end
            end
            % Setting on cid
            for idxCurrentCoord = newcid:obj.ncoords
                %cidをずらす．
                obj.cid(obj.dims_name(idxCurrentCoord)) = ...
                    obj.cid(obj.dims_name(idxCurrentCoord))+1;
            end
            obj.cid(dimName) = newcid;
            % Setting on coords
            obj.coords = [obj.coords(1:newcid-1), {coord(:)}, obj.coords(newcid:end)];
            % Setting on coords_name & dims_name
            obj.dims_name = [obj.dims_name(1:newcid-1), string(dimName), obj.dims_name(newcid:end)];
            obj.coords_name = [obj.coords_name(1:newcid-1), string(dimName), obj.coords_name(newcid:end)];
            % Setting on specified flag
            obj.coords_isSpecified = [obj.coords_isSpecified(1:newcid-1), isSpecified, obj.coords_isSpecified(newcid:end)];
            % Setting on dimension setting
            obj.coords_dims = [obj.coords_dims(:,1:newcid-1), false([obj.ncoords, 1]), obj.coords_dims(:,newcid:end)];
            obj.coords_dims = [obj.coords_dims(1:newcid-1,:); (1:newcid)==newcid; obj.coords_dims(newcid:end,:)];
            % increment of ndims and ncoords
            obj.ndims = obj.ndims +1;
            obj.ncoords = obj.ncoords +1;
            % どのような軸の追加であったとしても，既存のCoordinatesの中身には影響を与えない
        end
    end

    methods (Access=protected)
        function varargout = parenReference(obj, indexOp)
            assert(numel(obj.dims_name) == numel(indexOp(1).Indices) || numel(indexOp(1).Indices) == 1)
            % 各coordsのsliceの計算
            indices = indexOp(1).Indices;
            if(numel(indices)==1) %与えられたアドレスの値が一つだけ：2通り考えられる．
                if(obj.ndims > 1) %ndimsが2以上の場合は，(subsではなく）インデックス(ind)で与えたということ．
                    %インデックスをind2subを使って対応するsubsに変更する
                    indices = cell(1,obj.ndims);
                    [indices{:}] = ind2sub(obj.sz, indexOp(1).Indices{1});
                else %ndimsが1のときは，インデックスだが，縦ベクトルになるようにしなければならない
                    indices = cell(1,2); %縦ベクトルになるようにｓｕｂs指定できるようにする．
                    indices(1) = indexOp(1).Indices(1); %縦ベクトルにする
                    indices(2) = {1}; %縦ベクトルになる
                end
            end
            for idxCoords = 1:obj.ncoords
                subs = indices(obj.coords_dims(idxCoords,:));
                obj.coords{idxCoords} = obj.coords{idxCoords}(subs{:});
            end
            % データ配列のslice
%             obj.cont_array = obj.cont_array.(indexOp(1));
            obj.cont_array = obj.cont_array(indices{:});
            obj.sz = obj.csize();
            if isscalar(indexOp)
                varargout{1} = obj;
                return;
            end
            [varargout{1:nargout}] = obj.(indexOp(2:end));
        end

        function obj = parenAssign(obj,indexOp,varargin)
            % Ensure object instance is the first argument of call.
            if isempty(obj)
                obj = varargin{1};
            end
            if isscalar(indexOp)
                assert(nargin==3);
                rhs = varargin{1};
                obj.cont_array = xarray.DataArray.extendArray(obj.cont_array, indexOp.Indices);
                if(~(obj.ndims == 1 && size(obj.cont_array,2) == 1) && ...
                        obj.ndims < numel(size(obj.cont_array)))
                % dimensionが追加された場合：
                    for idxAdditionalDim = obj.ndims+1:numel(size(obj.cont_array))
                        obj.sz = obj.csize();
                        obj = obj.expand_single_dim('', idxAdditionalDim, 1:obj.sz(idxAdditionalDim), false);
                    end
                end
                obj.cont_array.(indexOp) = double(rhs);
                obj.sz = obj.csize();
                for idxCoord = 1:obj.ncoords
                    obj.coords{idxCoord} = xarray.DataArray.extendArray(obj.coords{idxCoord},num2cell(obj.sz(obj.coords_dims(idxCoord,:))));
                end
                %obj.CheckConsistency();
                return;
            end
            [obj.(indexOp(2:end))] = varargin{:};
        end

        function n = parenListLength(obj,indexOp,ctx)
            if numel(indexOp) <= 2
                n = 1;
                return;
            end
            containedObj = obj.(indexOp(1:2));
            n = listLength(containedObj,indexOp(3:end),ctx);
        end

        function obj = parenDelete(obj,indexOp)
            obj.cont_array.(indexOp) = [];
            % 各coordsのsliceの計算
            for idxCoords = 1:obj.ncoords
                subs = indexOp(1).Indices(obj.coords_dims(idxCoords,:));
                if(~all(strcmp(subs,':')))
                    obj.coords{idxCoords}(subs{:}) = [];
                end
            end
            obj.sz = obj.csize();
        end
        
    end
    
    % mathematical methods
    methods (Access=public)
        function out = value(obj)
            out = obj.cont_array;
        end
        
        function c = double(obj)
            c = obj.cont_array;
        end
        
        function r = mtimes(obj1,obj2)
            obj1mat = double(obj1);
            obj2mat = double(obj2);
            assert(numel(size(obj1mat))==2, 'xarray:invalidvalue',...
                'obj1 is not a matrix nor vector nor scalar');
            assert(numel(size(obj2mat))==2, 'xarray:invalidvalue',...
                'obj2 is not a matrix nor vector nor scalar');
            assert(isscalar(obj1mat) || isscalar(obj2mat) || size(obj1mat,2) == size(obj2mat,1),...
                'xarray:invalidvalue',...
                'mtimes(obj1,obj2) size(obj1,2) ~= size(obj2,1)');
            prodmat = obj1mat*obj2mat;
            if(isa(obj1, 'xarray.DataArray') && size(obj2mat,1) == size(obj2mat,2))
                r = obj1;
                r.cont_array = prodmat;
            elseif(isa(obj2, 'xarray.DataArray') && size(obj1mat,1) == size(obj1mat,2))
                r = obj2;
                r.cont_array = prodmat;
            else
                r = prodmatt;
            end
        end
        
        function r = times(obj1, obj2)
            obj1mat = double(obj1);
            obj2mat = double(obj2);
            if(isa(obj1,'xarray.DataArray'))
                r = obj1;
            else
                r = obj2;
            end
            r.cont_array = obj1mat.*obj2mat;
        end

        function out = sum(obj)
            out = sum(obj.cont_array,"all");
        end
        
        function out = cat(dim,varargin)
            numCatArrays = nargin-1;
            newArgs = cell(numCatArrays,1);
            for ix = 1:numCatArrays
                if isa(varargin{ix},'DataArray')
                    newArgs{ix} = varargin{ix}.ContainedArray;
                else
                    newArgs{ix} = varargin{ix};
                end
            end
            out = DataArray(cat(dim,newArgs{:}));
        end

        function varargout = size(obj,varargin)
            %このsize参照に限っては実際のcont_arrayのサイズを返す
            [varargout{1:nargout}] = size(obj.cont_array,varargin{:});
        end
        
        function str = char(obj)
            s = cell(1,3); ind = 1;
            s(ind) = {'<xarray.DataArray ('}; ind = ind+1;
            fspec = '%s:%2d';
            s(ind) = {sprintf(fspec, obj.dims_name(1), obj.sz(1))}; ind = ind+1;
            for idx_dim = 2:obj.ndims
                s(ind) = {sprintf([', ', fspec], obj.dims_name(idx_dim), obj.sz(idx_dim))};
                ind = ind+1;
            end
            s(ind) = {')>\n'}; ind=ind+1;
            
            s(ind) = {xarray.DataArray.charPages(obj.cont_array, obj.sz)}; ind=ind+1;
            
            if(any(obj.coords_isSpecified))
                s(ind) = {'Coordinates:\n'}; ind=ind+1;
                for idxCoord = 1:obj.ncoords
                    if(false == obj.coords_isSpecified(idxCoord))
                        continue;
                    end
                    % marker of dimensional coordinate
                    if(any(strcmp( obj.dims_name, obj.coords_name(idxCoord))))
                        s(ind)= {'  * '}; ind=ind+1;
                    else
                        s(ind) = {'    '}; ind=ind+1;
                    end
                    % name of coordinate
                    s(ind) = {sprintf('%-9s', obj.coords_name(idxCoord))}; ind = ind+1;
                    if(any(obj.coords_dims(idxCoord,:)))
                        s(ind) = {'('}; ind = ind+1;
                        s(ind) = {xarray.DataArray.charList('%s', obj.dims_name(obj.coords_dims(idxCoord,:)), ', ')}; ind = ind+1;
                        s(ind) = {')'}; ind = ind+1;
                    end
                    % class name of coordinate
                    s(ind) = {sprintf(' %s ', class(obj.coords{idxCoord}))}; ind = ind+1;
                    % list of coordinate values
                    s(ind) = {xarray.DataArray.charValInLine(obj.coords{idxCoord}, 60)}; ind = ind+1;
                    s(ind) = {'\n'}; ind = ind+1;
                end
            end

            if(~all(obj.coords_isSpecified(1:obj.ndims)))
                s(ind) = {'Dimensions without coordinates: '}; ind = ind+1;
                s(ind) = {xarray.DataArray.charList('%s', obj.dims_name(~obj.coords_isSpecified(1:obj.ndims)),', ')}; ind = ind+1;
                s(ind) = {'\n'}; ind = ind+1;
            end
            str = [s{1:ind-1}];
        end
        
        function disp(obj)
            fprintf(char(obj));
        end
    end

    methods(Access = private)
        function CheckConsistency(obj)%整合性のチェック
            %各coordに設定された要素の数が次元の数と整合しているかどうかのチェック
            for idxCoords = 1:obj.ncoords
                %coordinatesに登録されている配列のサイズ
                coordSize = size(obj.coords{idxCoords});
%                 coordNumel = prod(coordSize);
                %coordSizeの1の要素を除く
                coordSize = coordSize(coordSize~=1);
                %dimensionから計算される配列のサイズ
                dimSize = obj.sz(obj.coords_dims(idxCoords,:));
%                 dimNumel = prod(dimSize);
                %dimSizeの1の要素を除く
                dimSize = dimSize(dimSize~=1);
                assert(numel(coordSize) == numel(dimSize),...
                    "xarray:invalidcoorddims",...
                    "size('%s')=%s is different from %s",...
                    obj.coords_name(idxCoords), ...
                    xarray.DataArray.charList('%d', coordSize, 'x'),...
                    xarray.DataArray.charList('%d', dimSize, 'x')...
                    );
                assert(all(coordSize == dimSize),...
                    "xarray:invalidcoorddims",...
                    "size('%s')=%s is different from %s",...
                    obj.coords_name(idxCoords), ...
                    xarray.DataArray.charList('%d', coordSize, 'x'),...
                    xarray.DataArray.charList('%d', dimSize, 'x')...
                    )
            end
        end

        function sz = csize(obj)
            sz = [size(obj.cont_array), ones(1,obj.ndims-numel(size(obj.cont_array)))];
        end
    end

    methods (Static, Access=public)
        function obj = empty(varargin)
            obj = xarray.DataArray(double.empty(varargin{:}));
        end
    end


    methods(Static, Access = private)
        function str = charList(formatSpec, list, separator)
            str = string.empty;
            if(length(list) >= 1)
                str = sprintf( formatSpec, list(1));
                for idx = 2:length(list)
                    str = [str, sprintf( [separator, formatSpec] , list(idx))]; %#ok<AGROW> 
                end
            end
        end
        function str = charValInLine(vals,maxlength)
            s = cell(1,3); ind = 1;
            remainLength = maxlength - 3;
            printFormer = 1;
            printLatter = numel(vals);
            remainLength = remainLength - strlength(sprintf("%g%g", vals(printFormer), vals(printLatter)));
            while(printFormer < printLatter)
                if(remainLength > strlength(sprintf("%g",vals(printFormer+1))) + 1)
                    remainLength = remainLength - (strlength(sprintf("%g",vals(printFormer+1))) + 1);
                    printFormer = printFormer+1;
                else
                    break;
                end
                if(remainLength > strlength(sprintf("%g",vals(printLatter-1))) + 1)
                    remainLength = remainLength - (strlength(sprintf("%g",vals(printLatter-1))) + 1);
                    printLatter = printLatter-1;
                else
                    break;
                end
            end
            if(printFormer+1 < printLatter)
                s(ind) = {xarray.DataArray.charList('%g', vals(1:printFormer), ' ')};
                ind = ind+1;
                s(ind) = {' ... '}; ind = ind+1;
                s(ind) = {xarray.DataArray.charList('%g', vals(printLatter:end), ' ')};
                ind = ind+1;
            else
                s(ind) = {xarray.DataArray.charList('%g', vals(1:end), ' ')};
                ind = ind+1;
            end
            str = [s{1:ind-1}];
        end
        function str = charMatrix(vals, szVals)
            s = cell(1,3); ind = 1;
            if(szVals(1) > 5)
                s(ind) = {'\t['}; ind=ind+1;
                s(ind) = {xarray.DataArray.charValInLine(vals(1,:),80)}; ind = ind+1;
                s(ind) = {'\n\t '}; ind=ind+1;
                s(ind) = {xarray.DataArray.charValInLine(vals(2,:),80)}; ind = ind+1;
                s(ind) = {'\n\t '}; ind=ind+1;
                s(ind) = {'...\n\t '}; ind=ind+1;
                s(ind) = {xarray.DataArray.charValInLine(vals(end-1,:),80)}; ind=ind+1;
                s(ind) = {'\n\t '}; ind=ind+1;
                s(ind) = {xarray.DataArray.charValInLine(vals(end,:),80)}; ind=ind+1;
                s(ind) = {']\n'}; ind=ind+1;
            elseif(szVals(1) >= 1)
                s(ind) = {'\t['}; ind=ind+1;
                s(ind) = {xarray.DataArray.charValInLine(vals(1,:),80)};ind=ind+1;
                for i = 2:szVals(1)
                    s(ind) = {'\n\t '}; ind=ind+1;
                    s(ind) = {xarray.DataArray.charValInLine(vals(i,:),80)}; ind=ind+1;
                end
                s(ind) = {']\n'}; ind=ind+1;
            else
                s(ind) = {xarray.DataArray.charEmpty(vals)}; ind=ind+1;
            end
            str = [s{1:ind-1}];
        end
        function str = charPage(vals, szVals, pageIdx)
            s = cell(1,3); ind = 1;
            pageNums = szVals(3:end);
            pageSubs = cell(numel(pageNums),1);
            [pageSubs{:}] = ind2sub(pageNums, pageIdx);
            s(ind) = {'(:,:'}; ind=ind+1;
            for i=1:numel(pageNums)
                s(ind) = {sprintf(',%d',pageSubs{i})}; ind = ind+1;
            end
            s(ind) = {') =\n'}; ind=ind+1;
            s(ind) = {xarray.DataArray.charMatrix(vals(:,:, pageIdx), szVals)}; ind =ind+1;
            str = [s{1:ind-1}];
        end
        function str = charPages(vals, szVals)
            s = cell(1,3); ind = 1;
            if(numel(vals) == 0)
                str = xarray.DataArray.charEmpty(vals);
                return;
            end
            if(numel(szVals)<=2)
                str = xarray.DataArray.charMatrix(vals, szVals);
                return;
            end
            pageNums = szVals(3:end);
            numPages = prod(pageNums);
            if(numPages > 5)
                s(ind) = {xarray.DataArray.charPage(vals, szVals, 1)}; ind=ind+1;
                s(ind) = {xarray.DataArray.charPage(vals, szVals, 2)}; ind=ind+1;
                s(ind) = {'...\n'}; ind=ind+1;
                s(ind) = {xarray.DataArray.charPage(vals, szVals, numPages-1)}; ind=ind+1;
                s(ind) = {xarray.DataArray.charPage(vals, szVals, numPages)}; ind=ind+1;
            else
                for i = 1:numPages
                    s(ind) = {xarray.DataArray.charPage(vals, szVals, i)}; ind=ind+1;
                end
            end
            str = [s{1:ind-1}];
        end
        function str = charEmpty(vals)
            str = ['\tempty array:', xarray.DataArray.charList('%d',size(vals),'x'), ' ', class(vals), '\n'];
        end
        function arr = extendArray(arr, newIndices)
            currentSize = size(arr);
            newBoundary = cellfun(@max, newIndices);
            dim = max([numel(currentSize), numel(newBoundary)]);
            currentSize = [currentSize ones(1,dim-numel(currentSize))];
            newBoundary = [newBoundary ones(1,dim-numel(newBoundary))];
            dims = 1:dim;
            for idxDim = dims(currentSize < newBoundary)
                concatsize = size(arr);
                concatsize(idxDim) = newBoundary(idxDim)-currentSize(idxDim);
                arr = cat(idxDim, arr, NaN(concatsize));
            end
        end
    end
end

