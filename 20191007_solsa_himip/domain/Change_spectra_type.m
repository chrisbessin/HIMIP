classdef Change_spectra_type
   % Modification of the reflectance
   %
   % wavelengths : list of double
   %           wavelengths used in the spectra
   % raw_spectra : list of double
   %           spectra to be modified, it has to be a relative reflectance
   % 
   properties
      wavelengths
      raw_spectra
      spectrum_type
   end

   methods
        function obj = Change_spectra_type(wavelengths,raw_spectra)
           % Creation of the object change spectra
           %
           %
            obj.wavelengths = wavelengths;
            obj.raw_spectra = raw_spectra;
        end
        function spectra = apply(obj,spectrum_type)
            % Modification the reflectance by apply a log10 or removing
            % the hull 
            %
            % Parameters
            % ----------
            % spectrum_type : int
            %           output type of spectrum
            %
            % Output
            % ------
            % spectra : list of double
            %           output spectra after modification
            %
            
            if spectrum_type == 0           
                % No modification
                spectra = obj.raw_spectra;
            
            elseif spectrum_type==1       
                % Log10 data
                spectra = log10(obj.raw_spectra + 1.0) + eps;
            
            elseif spectrum_type==2
                % Noise removal
                spectra_temp = obj.raw_spectra;
                spectra_temp([false(size(spectra_temp, 1), 1),...
                              diff(spectra_temp, 2, 2) > 0.1,...
                              false(size(spectra_temp, 1), 1)]) = NaN;
                spectra_temp = fillmissing(spectra_temp,...
                                           'linear', 2);
                spectra_temp = loop_over_spectra(spectra_temp,...
                                                 'sgolayfilt', 1, 13);
                spectra = spectra_temp;
            
            elseif spectrum_type==3        
                % Continuum removal
                spectra_temp = obj.raw_spectra;
                spectra_temp([false(size(spectra_temp, 1), 1),...
                              diff(spectra_temp,2,2) > 0.1,...
                              false(size(spectra_temp, 1), 1)]) = NaN;
                spectra_temp = fillmissing(spectra_temp,...
                                           'linear',2);
                spectra_temp = loop_over_spectra(spectra_temp,...
                                                 'ContinuumRemoval',...
                                                 obj.wavelengths,...
                                                 false);
                spectra = spectra_temp;
            
            elseif spectrum_type==4
                % hull (trend)
                spectra_temp = obj.raw_spectra;
                spectra_temp([false(size(spectra_temp, 1), 1),...
                              diff(spectra_temp, 2, 2) > 0.1,...
                              false(size(spectra_temp, 1), 1)]) = NaN;
                spectra_temp = fillmissing(spectra_temp,...
                                           'linear', 2);
                spectra_temp = loop_over_spectra(spectra_temp,...
                                                 'ContinuumRemoval',...
                                                 obj.wavelengths,...
                                                 true);
                spectra = spectra_temp;
            
            elseif spectrum_type == 5
                % Noise removal + Continuum removal
                spectra_temp = obj.raw_spectra;
                spectra_temp([false(size(spectra_temp, 1), 1),...
                              diff(spectra_temp, 2, 2) > 0.1,...
                              false(size(spectra_temp, 1), 1)]) = NaN;
                spectra_temp = fillmissing(spectra_temp,'linear', 2);
                spectra_temp = loop_over_spectra(spectra_temp,...
                                                 'sgolayfilt', 1, 13);
                spectra_temp = loop_over_spectra(spectra_temp,...
                                                 'ContinuumRemoval',...
                                                 obj.wavelengths,...
                                                 false);
                spectra = spectra_temp;
            
            elseif spectrum_type == 6
                % Continuum removal + Noise removal
                spectra_temp = obj.raw_spectra;
                spectra_temp([false(size(spectra_temp, 1), 1),...
                              diff(spectra_temp, 2, 2) > 0.1,...
                              false(size(spectra_temp, 1), 1)]) = NaN;
                spectra_temp = fillmissing(spectra_temp,...
                                           'linear', 2);
                spectra_temp = loop_over_spectra(spectra_temp,...
                                                 'ContinuumRemoval',...
                                                 obj.wavelengths,...
                                                 false);
                spectra_temp = loop_over_spectra(spectra_temp,...
                                                 'sgolayfilt', 1, 13);
                spectra = spectra_temp;
            
            elseif spectrum_type == 7
                % Derivative order 1
                spectra_temp = obj.raw_spectra;
                spectra_temp([false(size(spectra_temp, 1), 1),...
                              diff(spectra_temp, 2, 2) > 0.1,...
                              false(size(spectra_temp, 1), 1)]) = NaN;
                spectra_temp = fillmissing(spectra_temp,...
                                           'linear',2);
                spectra_temp = loop_over_spectra(spectra_temp,...
                                                 'sgolayfilt', 1, 61);
                spectra_temp = loop_over_spectra(spectra_temp,...
                                                 'derivative', 1);
                spectra_temp = loop_over_spectra(spectra_temp,...
                                                 'sgolayfilt', 1, 15);
                spectra_temp = loop_over_spectra(spectra_temp,...
                                                 'normalise', 0.5, []); %[-0.002,0.002]
                %spectra_temp = spectra_temp > 0.5;
                spectra = spectra_temp;
                
            elseif spectrum_type==8
                % Derivative order 2
                spectra_temp = obj.raw_spectra;
                spectra_temp([false(size(spectra_temp, 1), 1),...
                              diff(spectra_temp, 2, 2) > 0.1,...
                              false(size(spectra_temp, 1), 1)]) = NaN;
                spectra_temp = fillmissing(spectra_temp,...
                                           'linear',2);
                spectra_temp = loop_over_spectra(spectra_temp,...
                                                 'sgolayfilt', 1, 41);
                spectra_temp_1 = loop_over_spectra(spectra_temp,...
                                                   'derivative', 1);
                spectra_temp = loop_over_spectra(spectra_temp_1,...
                                                 'sgolayfilt', 1, 41);
                spectra_temp = loop_over_spectra(spectra_temp,...
                                                 'derivative', 1);
                spectra_temp = loop_over_spectra(spectra_temp,...
                                                 'sgolayfilt', 1, 15);
                spectra_temp = spectra_temp.*spectra_temp_1;
                spectra_temp = loop_over_spectra(spectra_temp,...
                                                 'sgolayfilt', 1, 15);
                spectra_temp = loop_over_spectra(spectra_temp,...
                                                 'normalise', 0.5, []); %[-0.00000005,0.00000005]
                %spectra_temp = spectra_temp > 0.75;
                spectra = spectra_temp;
            
            elseif spectrum_type==9
                % Bollinger bands
                spectra_temp = obj.raw_spectra;
                spectra_temp([false(size(spectra_temp, 1), 1),...
                              diff(spectra_temp, 2, 2) > 0.1,...
                              false(size(spectra_temp, 1), 1)]) = NaN;
                spectra_temp = fillmissing(spectra_temp,...
                                           'linear', 2);
                spectra_temp = loop_over_spectra(spectra_temp,...
                                                 'ContinuumRemoval',...
                                                 obj.wavelengths,...
                                                 false);
                spectra_temp = loop_over_spectra(spectra_temp,...
                                                 'bollinger', 15, 0.5);
                spectra = spectra_temp;
            
            elseif spectrum_type==10
                spectra_temp = obj.raw_spectra;
                spectra_temp = loop_over_spectra(spectra_temp,...
                                                 'sgolayfilt',1,3);
                max_values = max(spectra_temp, [], 2);
                spectra_temp = spectra_temp./max_values;
                spectra = spectra_temp;
            
            else
                error('Invalid spectrum type')
            end
        end
   end
end

function spectra = loop_over_spectra(input_spectra,funct_app,option1,option2)
    % This function performs hull corrections for the input reflectance, 
    % i.e., removes the background
    %
    % Parameters
    % ----------
    % raw_spectra: array of double (1 x n)
    %       raw_spectra data
    % funct_app: string
    %       transformation applied to the spectra
    % option1: 
    %       1st option of the transformation
    % option2: 
    %       2nd option of the transformation
    % 
    % Output
    % ------
    % spectra: array of double (1 x n)
    %       the continuum removal of the input spectrum
    %
    
    % Check if the dimensions of the spectra, if it's 3D, then it is
    % reshaped to a 2D matric to allow processing
    if length(size(input_spectra))==3
        [NB_samples, NB_lines, NB_wls] = size(input_spectra);
        NB_spectra = NB_lines*NB_samples;
        raw_spectra = reshape(input_spectra,[NB_spectra,NB_wls]);
    else
        [NB_spectra,NB_wls] = size(input_spectra);
        raw_spectra = input_spectra;
    end
    
    % Output variable is set
    spec_temp = zeros([NB_spectra, NB_wls]);

    % Selection of the lines 
    ind_not_nan_data = ~isnan(input_spectra);
    ind_sel_lin = max(ind_not_nan_data,[],2)==1;
    index_lin = 1:NB_spectra;
    index_not_nan = index_lin(ind_sel_lin);

    for j = 1:length(index_not_nan)
        i = index_not_nan(j);

        % Application of the transformation
        if strcmp(funct_app,'sgolayfilt')
            spec_temp(i,:) = sgolayfilt(raw_spectra(i,:),option1,option2);

        elseif strcmp(funct_app,'derivative')
            temp_1 = squeeze(diff(raw_spectra(i,:),option1));
            even_part = floor(option1/2);   % Add zeros symmetrically on both ends of the spectrum
            temp_2 = [zeros(1,even_part), temp_1, zeros(1,even_part + mod(option1,2))];
            spec_temp(i,:) = temp_2;%.*squeeze(raw_spectra(i,j,:))';

        elseif strcmp(funct_app,'ContinuumRemoval')
            spec_temp(i,:) = ContinuumRemoval(raw_spectra(i,:),option1,option2);

        elseif strcmp(funct_app,'bollinger')
            spec_temp(i,:) = bollinger_bands(raw_spectra(i,:),option1,option2);

        elseif strcmp(funct_app,'normalise')
            spec_temp(i,:) = norm_center(raw_spectra(i,:),option1,option2);
        else
           error('Transformation not known...')

        end
    end
    spectra = spec_temp;
    
    % If input is dimension 3, then it is reshaped back to its original
    % format
    if length(size(input_spectra))==3
        spectra = reshape(spectra,[NB_samples, NB_lines, NB_wls]);
    end
    
end

function spectra = norm_center(raw_spectra,centered_val,norm)
    % This function performs hull corrections for the input reflectance, 
    % i.e., removes the background
    %
    % Parameters
    % ----------
    % raw_spectra: array of double (1 x n)
    %       the relative reflectance of hyperspectral data
    % centered_val: double
    %       value around which the spectra will be centered
    % max_val: double
    %       absolute maximum value, the spectra will be divide by that
    %       value
    %     
    % Output
    % ------
    % spectra: array of double (1 x n)
    %       normalised and centered spectra
    %
    
    if isempty(norm)
        norm = max(abs(raw_spectra));
    elseif length(norm)==2
        norm = max(norm) - min(norm);
    end
    
    spectra =  raw_spectra/norm/2 + centered_val;
    
end

function spectra = ContinuumRemoval(raw_spectra, wavelengths,hull_trend)
    % This function performs hull corrections for the input reflectance, 
    % i.e., removes the background
    %
    % Parameters
    % ----------
    % wavelengths: array of double (1 x n)
    %       wavelengths of the reflectance
    % raw_spectra: array of double (1 x n)
    %       the relative reflectance of hyperspectral data
    % hull_trend: boolean
    %       option to save the hull function
    %     
    % Output
    % ------
    % spectra: array of double (1 x n)
    %       the continuum removal of the input spectrum
    %
    % Author: Thanh Bui (thanh.bui@erametgroup.com)
    
    % NaN management
    idx_nan = isnan(raw_spectra); 
    idx_no_nan = ~idx_nan;
    raw_spec_no_nan = raw_spectra(idx_no_nan);
    wls_no_nan = wavelengths(idx_no_nan);
    spectra(idx_nan) = NaN(1,sum(idx_nan));
    
    % Compute hull curve and hull removal 
    if mean(diff(diff(raw_spec_no_nan)))~=0 % if the input is not a line
        [hull,~] = convhull(wls_no_nan, raw_spec_no_nan);
        start_idx = find(hull == length(raw_spec_no_nan));
        up_hull = hull(end:-1:start_idx);
        refl_hull = interp1(wls_no_nan(up_hull), raw_spec_no_nan(up_hull), wls_no_nan);
    else
        refl_hull = raw_spec_no_nan;
    end

    % Outputs
    if hull_trend
        spectra(idx_no_nan) = refl_hull;
    else
        refl_cr = raw_spec_no_nan./(refl_hull + eps);
        spectra(idx_no_nan) = refl_cr;
    end
end

function spectra = bollinger_bands(raw_spectra, size_mov_window, Nb_std)
    %
    %
    %
    raw_spectra = squeeze(raw_spectra);
    raw_spectra = raw_spectra';
    
    mov_ave = movmean(raw_spectra,size_mov_window);
    mov_std = movstd(raw_spectra,int16(size_mov_window*1.));
    
    %upper = mov_ave + Nb_std*mov_std;
    lower = mov_ave - Nb_std*mov_std;
    
    spectra = (raw_spectra < lower);
    
    %{
    figure(2);
    hold on
    plot(raw_spectra)
    plot(mov_ave)
    plot(lower)
    plot(upper)
    %}
end