function [V, options] = christiaans_v(csv_path, options)

tic
k = 1;
fid = fopen(csv_path, 'r');
clear('X')
while ~feof(fid)
    line    = fgetl(fid);
    y       = textscan(line,'%s', 'Delimiter', ',', 'TreatAsEmpty',  ' ');
    varnames(k)  = string(y{1}{1}); %#ok<SAGROW>
    X(k).Val   = str2double(y{:});
    k       = k+1;
end
fclose(fid);
toc
V                           = assignvarnames();

for i = 1:length(V)
    j = find(strcmp(varnames,V(i).Name));
    if isempty(j)
        if i==2
            fprintf(1,'%s %s %s \n','warning: input "', V(i).Name, '" not provided in input data...');
            fprintf(1,'%s %s %s\n', 'I will use 0.25*Cab instead');
            options.Cca_function_of_Cab = 1;
        else
            if ~(options.simulation==1) && (i==30 || i==32)
                fprintf(1,'%s %s %s \n','warning: input "', V(i).Name, '" not provided in input data...');
                fprintf(1,'%s %s %s\n', 'I will use the MODTRAN spectrum as it is');
            else
                if (options.simulation == 1 || (~options.simulation && (i<46 || i>50 ))) && i<65
                    fprintf(1,'%s %s %s \n','warning: input "', V(i).Name, '" not provided in input data');
                    if (options.simulation ==1 && (i==1 ||i==9||i==22||i==23||i==54 || (i>29 && i<37)))
                        fprintf(1,'%s %s %s\n', 'I will look for the values in Dataset Directory "',F(5).FileName,'"');
                    else
                        if (i== 24 || i==25)
                            fprintf(1,'%s %s %s\n', 'will estimate it from LAI, CR, CD1, Psicor, and CSSOIL');
                            options.calc_zo = 1;
                        else
                            if (i>38 && i<44)
                                fprintf(1,'%s %s %s\n', 'will use the provided zo and d');
                                options.calc_zo = 0;
                            else
                                if ~((options.simulation ==1 && (i==30 ||i==32)))
                                    fprintf(1,'%s \n', 'this input is required: SCOPE ends');
                                    return
                                elseif (options.simulation ==1 && (i==30 ||i==32))
                                    fprintf(1,'%s %s %s\n', '... no problem, I will find it in Dataset Directory "',F(5).FileName, '"');
                                end
                            end
                        end
                    end
                end
            end
        end
    else
        k = find(~isnan(X(j).Val));
        if ~isempty(k)
            V(i).Val = X(j).Val(k);
        else
            V(i).Val            = -999;
        end
    end
end
end