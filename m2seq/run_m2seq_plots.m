function run_m2seq_plots( r_dms, r_nomod, name, output_folder, iterations, offset, shapeknots, fiveprime_trim, threeprime_trim, new_struct);
% Adapted from Ved Topkar

    if ~exist(output_folder, 'dir')
	    mkdir(output_folder)
    end

    % Load in the rdat files
    rna_dms = read_rdat_file(r_dms);
    rna_nomod = read_rdat_file(r_nomod);

    % Extract the reference sequence from the DMS rdat
    sequence = rna_dms.sequence;
    sequence_length = length(sequence);
    trimmed_sequence = sequence(fiveprime_trim:sequence_length-threeprime_trim);

    % Compute the Zscore matrix
    Z = output_Zscore_from_rdat( [], r_dms, r_nomod, [], 1, 1 );
    % Subset Zscores to trim ends
    Z = Z(fiveprime_trim:sequence_length-threeprime_trim, fiveprime_trim:sequence_length-threeprime_trim)
    show_2d_das(-Z, name, offset, output_folder);

    % Get 1D dms profile (again, from dms and nomod rdats)
    one_dimensional_dms = get_DMS_profile(r_dms, r_nomod)
    one_dimentional_dms = one_dimensional_dms(fiveprime_trim:sequence_length-threeprime_trim)
    one_dimentional_dms_normalized = DMS_normalize(one_dimensional_dms, sequence)
    one_dimentional_dms_normalized = one_dimentional_dms_normalized(fiveprime_trim:sequence_length-threeprime_trim)

    % Save 1D data to CSV
    csvwrite('1D_DMS_normalized.csv', one_dimentional_dms_normalized)

    % Predict helices with M2Net
    clf
    close
    m2net_figure = figure
    [Z_m2net, structure_m2net] = m2net( Z,  trimmed_sequence, '', 0);
    saveas(gcf, fullfile(output_folder, [name, '_', 'M2net_plot', '.png']))
    clf
    close
    show_2dmap(Z_m2net, '', 0);
     
    % Run bootstrapped RNA_structure with X iterations
    bootstrapped_figure = figure
    [final_structure, bpp, one_dimensional_filtered] = rna_structure(trimmed_sequence, '', '', '', Z, iterations, shapeknots, one_dimentional_dms_normalized);
    clf
    close

    % Generate figure for bootstraped pairing probabilities
    bpp_figure = figure
    show_2dmap(bpp, '', offset, 1, 0);
    pbaspect([1 1 1]); % Square aspect ratio
    set(gca, 'FontSize', 14);
    xlabel('Sequence Position', 'FontSize', 20, 'FontWeight', 'bold');
    ylabel('Sequence Position', 'FontSize', 20, 'FontWeight', 'bold');
    colorbar;
    title('Bootstrapped Base Pairing Probabilities');
    saveas(gcf, fullfile(output_folder, [name, '_', 'bpp', '.png']));
    savefig(fullfile(output_folder, [name, '_bpp','.fig']));
    clf
    close

    % Generate VARNA fig
    % We pass in bpp to display the bootstrapping-derived helix probabilities
    figure
    output_varna(fullfile(output_folder, [name, '_VARNA_pred' ,'.png']), trimmed_sequence, final_structure, '', '', 0, '', [], one_dimentional_dms_normalized, bpp, '')
    output_varna(fullfile(output_folder, [name, '_VARNA_pred']), trimmed_sequence, final_structure, '', '', 0, '', [], one_dimentional_dms_normalized, bpp, '')

    if ~strcmp(new_struct, '');
        figure
        output_varna(fullfile(output_folder, [name, '_VARNA' ,'.png']), trimmed_sequence, new_struct, '', '', 0, '', [], one_dimentional_dms_normalized, bpp, '')
        output_varna(fullfile(output_folder, [name, '_VARNA']), trimmed_sequence, new_struct, '', '', 0, '', [], one_dimentional_dms_normalized, bpp, '')
    end
end


function show_2d_das(bpp, name, offset, output_folder)
    show_2dmap(bpp, '', offset, 1, 0);
    pbaspect([1 1 1]); % Square aspect ratio
    set(gca, 'FontSize', 14);
        xlabel('Sequence Position', 'FontSize', 20, 'FontWeight', 'bold');
        ylabel('Mutation Position', 'FontSize', 20, 'FontWeight', 'bold');
    saveas(gcf, fullfile(output_folder, [name, '_Z.pdf']))
end