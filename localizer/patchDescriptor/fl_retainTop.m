function binMaskTop = retainTop( response, binMask, NUM)
binMaskTop = zeros( size(binMask,1), size(binMask,2) );
[y x]= find(binMask ~= 0 );
allVals = response( binMask ~= 0 );
[sortedRes sortedIds] = sort( allVals, 'descend' );
NUM = min( NUM, numel(y) );
for iter = 1:NUM
    binMaskTop( y( sortedIds(iter) ), x( sortedIds(iter) ) ) = 1;
end