% For each monkey and segment, compute the row-wise mean across and store in a separate output file
    % Access the 'int' matrix for current data_monks{1, idx}
    int_matrix = data_monks{1, idx}.int; % n x 4 x 5 matrix
    name = data_monks{1, idx}.monkey; % Monkey name as a string
    num_rows = 4158;
    
    
    % Loop through each of the 4 "segments" (columns in the middle dimension)
    for segment = 1:4
        % Extract the slice corresponding to the current segment
        segment_matrix1 = int_matrix(:, segment, :); % n x 1 x 5 matrix
        
        % Squeeze to remove singleton dimensions (n x 5 matrix)
        segment_matrix2 = squeeze(segment_matrix1);
        segment_matrix = segment_matrix2(crd_acc);
       
        % Use the first column of the averaged data as the replacement for column 4 in data_acc
        row_avg = mean(segment_matrix, 2);
        data_acc(1:num_rows, 4) = row_avg(1:num_rows);
        
        % Save the modified data_acc for this segment to a separate ASCII file
        filename = sprintf('dlPFC_%s_segment%d', name, segment);
        save(filename, 'data_acc');
    end
end
