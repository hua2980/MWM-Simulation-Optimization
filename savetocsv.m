dims = size(PMs);

outf = [];
if (length(dims) == 3)
    for i = 1:dims(1)
        for j = 1:dims(2)
            for k = 1:dims(3)
                outf = [outf; i j k PMs(i,j,k)];
            end
        end
    end
elseif (length(dims) == 4)     
    for i = 1:dims(1)
        for j = 1:dims(2)
            for k = 1:dims(3)
                for l = 1:dims(4)
                    outf = [outf; i j k l PMs(i,j,k,l)];% Pars(j,1:4)];
                end
            end
        end
    end
end
csvwrite('etdecay95.csv',outf);