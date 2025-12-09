% 生成分组情况
function groups = generateGroups(n, k)
    groups = cell(stirling(n, k), 1);
    indices = 1:n;
    currentGroup = 1;
    
    generateGroup(indices, n, k, currentGroup);
    
    function generateGroup(indices, n, k, groupNum)
        if k == 0
            groups{groupNum} = indices(1:(n-k));
            return
        end
        
        for i = 1:(n-k+1)
            subGroup = indices(i);
            generateGroup(indices((i+1):end), n-i, k-1, groupNum+1);
            if groupNum <= size(groups, 1)
                groups{groupNum} = [subGroup, groups{groupNum}];
            end
        end
    end
end
