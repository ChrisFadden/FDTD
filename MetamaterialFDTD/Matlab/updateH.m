function [field] = updateH(field,grid,material)
    %Update Equations for FDTD
    
    %Free Space
    for i = 1:material(1).start - 1
        field.Hy(i) = field.Chyh*field.Hy(i) + field.Chye*(field.Ez(i) - field.Ez(i+1));
    end
    
     %Material
     Chyh = (1 - material(1).magneticLoss) / (1 + material(1).magneticLoss);
     Chye = (field.Chye / material(1).mur) / (1 + material(1).magneticLoss);
     for i = material(1).start:material(1).end       
        field.Hy(i) = Chyh*field.Hy(i) + Chye*(field.Ez(i) - field.Ez(i+1));
     end
    
    %Materials
    for m = 2:numel(material)
        %Free Space
        for i = material(m-1).end+1:material(m).start-1
            field.Hy(i) = field.Chyh*field.Hy(i) + field.Chye*(field.Ez(i) - field.Ez(i+1));
        end
        
        %Material
        Chyh = (1 - material(m).magneticLoss) / (1 + material(m).magneticLoss);
        Chye = (field.Chye / material(m).mur) / (1 + material(m).magneticLoss);
        for i = material(m).start:material(m).end       
            field.Hy(i) = Chyh*field.Hy(i) + Chye*(field.Ez(i) - field.Ez(i+1));
        end
    end  
    
    %Free Space
    for i = material(numel(material)).end+1:grid.sizeX-1
        field.Hy(i) = field.Chyh*field.Hy(i) + field.Chye*(field.Ez(i) - field.Ez(i+1));
    end
    
end


