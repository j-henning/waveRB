

height = 4;
tree = cell(2^height-1,1);
for i = 1:2^height - 1
    tree{i} = struct('muMin', NaN, 'muMax', NaN, 'U', NaN, 'Y', 'NaN', ...
        'S', 'NaN', 'snapshot', NaN);
end




for i=1:5
    mu = rand(1,1);
    index = 1; % start with the root
    while true
        if isnan(tree{index}.mu)
            tree{index}.mu = mu;
            break;   
        elseif mu < tree{index}.mu
            index = 2 * index; % Left child
        elseif mu > tree{index}.mu
            index = 2 * index + 1; % Right child
        else
            error('Parameter is already stored in the binary tree')
        end
        
        if index > 2^height-1
            error('Tree is not sufficiently large')
        end
    end
    
end



% Printing
linebreaks = 2.^(1:height)-1;
spacewidth = 16;
for i =1:2^height-1
    for j=1:spacewidth
        fprintf(' ')
    end
    fprintf('%.2f\t', tree{i}.mu)
    
    if ismember(i, linebreaks)
        fprintf('\n')
        spacewidth = spacewidth / 2;
    end
        
end

