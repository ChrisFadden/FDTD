function [field] = updateE(field,grid,material) 
    %Update Equations for FDTD 
    
    %Free Space
    for i = 2:material(1).start-1
       field.Ez(i) = field.Ceze*field.Ez(i) + field.Cezh*(field.Hy(i-1) - field.Hy(i)); 
    end
    
    %Material
    Ceze = (1 - material(1).electricLoss) / (1 + material(1).electricLoss);
    Cezh = (field.Cezh / material(1).epsr) / (1 + material(1).electricLoss);
    for i = material(1).start:material(1).end
       field.Ez(i) = Ceze*field.Ez(i) + Cezh*(field.Hy(i-1) - field.Hy(i)); 
    end
    
    for m = 2:numel(material)    
        %Free Space
        for i = material(m-1).end+1:material(m).start-1
            field.Ez(i) = field.Ceze*field.Ez(i) + field.Cezh*(field.Hy(i-1) - field.Hy(i)); 
        end
        
        %Material
        Ceze = (1 - material(m).electricLoss) / (1 + material(m).electricLoss);
        Cezh = (field.Cezh / material(m).epsr) / (1 + material(m).electricLoss);
        for i = material(m).start:material(m).end
            field.Ez(i) = Ceze*field.Ez(i) + Cezh*(field.Hy(i-1) - field.Hy(i)); 
        end
    end
     
    %Free Space
    for i = material(numel(material)).end+1:grid.sizeX
       field.Ez(i) = field.Ceze*field.Ez(i) + field.Cezh*(field.Hy(i-1) - field.Hy(i)); 
    end
    
end

