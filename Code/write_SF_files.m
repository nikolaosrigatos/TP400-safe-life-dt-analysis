function write_SF_files(TH_profile, toggle, num_permutations, folder_path)

    TH_profile_NASGRO = TH_profile; 
    TH_profile_NASGRO(:, 2) = TH_profile_NASGRO(:, 2) / 3;
    TH_profile_NASGRO(:, 3) = TH_profile_NASGRO(:, 3) / 3;

    switch toggle
        case 'yes'
            for i = 1:num_permutations
                % Construct the full file path
                new_filename = [folder_path, sprintf('SF_permuted_%d.txt', i)];
                
                % Check if the file exists before attempting to delete it
                if exist(new_filename, 'file')
                    delete(new_filename);
                end

                permuted_data = TH_profile_NASGRO(randperm(size(TH_profile_NASGRO, 1)), :);
                permuted_data = permuted_data(permuted_data(:, 1) ~= 0, :);

                fid = fopen(new_filename, 'w');
                fprintf(fid, '%d\t%.8f\t%.8f\n', permuted_data');
                fclose(fid);

            end
            
            fprintf('SF files were generated\n');
            fprintf('Fill the t_NASGRO_profile array\n');
            keyboard;

        case 'no'
            fprintf('SF files will not be generated');
    end

end