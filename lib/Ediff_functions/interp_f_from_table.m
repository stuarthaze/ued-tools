% function to return |f| and phase(f) for given s values, atomic
% numbers Z and energy by interpolating based on s and energy
% Input- Z, s: 1D array containing atomic numbers and s values resp.,
%        E: Energy in keV, pathname: path for folder with the scattering factors
% Output- |f| and phase as 2D arrays with the rows representing the s values 
%         and columns representing the elements correspoding to the Z array

function [mod_f, phase_f] = interp_f_from_table (Z, s, E, pathname) 
    % converting s and Z into column vectors
    if (~iscolumn(s))
        s = s.' ;
    end

    if (~iscolumn(Z))
        Z = Z.' ;
    end

    files = dir([pathname,'\*.txt']);  
    mod_f = zeros(length(s), length(Z));
    phase_f = zeros(length(s), length(Z));
    Energy_list = zeros(size(files,1), 1);
    searchPattern = '\d*(?=keV)';
    
    % Creating list of energies available in the scattering factors folder
    % and sorting it in ascending order
    for i = 1: size(files,1)
        energy = regexp(files(i).name, searchPattern, 'match');
        Energy_list(i) = str2num(energy{1});
    end

    Energy_list = sort(Energy_list);
    
    % Directly getting the scattering factors from the table if input
    % energy matches one of the energies in the folder
    for i = 1: length(Energy_list)
        if (E == Energy_list(i))
            fname= [pathname, '\SCATFAC_', num2str(E), 'keV.txt'];
            inputdata = load(fname);
            old_s = inputdata(:,1);
            for j = 1:length(Z)
                % getting the values baesd on atomic number
                old_mod_f = inputdata(:, Z(j)*2);
                old_phase_f = inputdata(:, Z(j)*2 + 1);
                % interpolating to get |f| and phase values for the input s
                new_mod_f = interp1 (old_s, old_mod_f, s, 'spline');
                new_phase_f = interp1 (old_s, old_phase_f, s, 'spline');
                % storing the values for each atom in different columns
                % based on the input Z vector
                mod_f(:, j) =  new_mod_f;
                phase_f (:,j) = new_phase_f;
            end
        end
    end
    
    % If input energy does not exactly match to any energies in the folder, 
    % interpolating linearly between the two energies it lies between
    for i = 1: length(Energy_list)-1
        % selecting the two energies it lies between
        if ( (E > Energy_list(i)) && (E < Energy_list(i+1)) )
            two_energies = [Energy_list(i) Energy_list(i+1)];
            fname1= [pathname, '\SCATFAC_', num2str(Energy_list(i)), 'keV.txt'];
            fname2= [pathname, '\SCATFAC_', num2str(Energy_list(i+1)), 'keV.txt'];
            inputdata1 = load(fname1);
            inputdata2 = load(fname2);
            old_s1 = inputdata1(:,1);
            old_s2 = inputdata2(:,1);
            for j = 1:length(Z)
                old_mod_f1 = inputdata1(:, Z(j)*2);
                old_phase_f1 = inputdata1(:, Z(j)*2 + 1);

                old_mod_f2 = inputdata2(:, Z(j)*2);
                old_phase_f2 = inputdata2(:, Z(j)*2 + 1);
                
                % First interpolating using spline based on input s 
                new_mod_f1 = interp1 (old_s1, old_mod_f1, s, 'spline');
                new_phase_f1 = interp1 (old_s1, old_phase_f1, s, 'spline');

                new_mod_f2 = interp1 (old_s2, old_mod_f2, s, 'spline');
                new_phase_f2 = interp1 (old_s2, old_phase_f2, s, 'spline');

                new_mod_f = [new_mod_f1 new_mod_f2];
                new_phase_f = [new_phase_f1 new_phase_f2];
                
                % Then interpolating based on input energy
                for k = 1: length(s)
                    mod_f(k, j) =  interp1(two_energies, new_mod_f(k,:), E);
                    phase_f (k,j) = interp1(two_energies, new_phase_f(k,:), E);
                end
            end
        end
    end
    
end
