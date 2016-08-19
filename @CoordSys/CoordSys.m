classdef CoordSys < handle
    properties
        num_cells
        arrays_sizes
        origin_cells=1
        
        lattice_const
        
        coord_stps
        coord_limits=[0,0];
        coord_sizes
        origin_inds
    end
    
    methods
        
        function obj = CoordSys(num_cells,arrays_sizes,units)
            
            if strcmp(units,'au')
                obj.lattice_const=0.5431*1e-9/MyConst.ab;
            else
                obj.lattice_const=0.5431*1e-9;
            end;
            
            obj.num_cells = num_cells;
            obj.arrays_sizes = arrays_sizes;            
            obj.coord_sizes = obj.num_cells*obj.arrays_sizes+1;
            
            dist=obj.num_cells*obj.lattice_const;
            obj.coord_stps=dist/(obj.coord_sizes-1);
            obj.origin_inds=obj.arrays_sizes*(obj.origin_cells-1)+1;
            obj.coord_limits(1)=-obj.origin_inds*obj.coord_stps;
            obj.coord_limits(2)=-obj.origin_inds*obj.coord_stps+obj.num_cells*obj.lattice_const;
            
        end
        
        function set_origin_cells(obj,origin_cells)
            obj.origin_cells = origin_cells;
            obj.origin_inds=obj.arrays_sizes*(obj.origin_cells-1)+1;
            obj.coord_limits(1)=-obj.origin_inds(1)*obj.coord_stps+obj.coord_stps;
            obj.coord_limits(2)=-obj.origin_inds(1)*obj.coord_stps+obj.num_cells*obj.lattice_const;
        end
        
%         [kx ,ky ,kz] = global_basis_rec(obj,j1,j2,j3)
%         [x ,y ,z] = global_basis(obj,j1,j2,j3)
%         
%         function [x ,y ,z] = global_basis_array(obj,X,Y,Z)
%             [x ,y ,z] = arrayfun(@obj.global_basis,X,Y,Z);
%         end;
%         
%         function [x ,y ,z] = global_basis_rec_array(obj,X,Y,Z)
%             [x ,y ,z] = arrayfun(@obj.global_basis_rec,X,Y,Z);
%         end;
%         
%         function [x ,y ,z] = global_coords_gen(obj)
%             x=1:obj.coord_sizes;
%             [X,Y,Z]=meshgrid(x,x,x);
%             [x ,y ,z] = arrayfun(@obj.global_basis,X,Y,Z);
%         end;
%         
%         function [x ,y ,z] = global_coords_rec_gen(obj)
%             x=1:obj.coord_sizes;
%             [X,Y,Z]=meshgrid(x,x,x);
%             [x ,y ,z] = arrayfun(@obj.global_basis_rec,X,Y,Z);
%         end;
        
        function x = x(obj)
            x = obj.coord_limits(1):obj.coord_stps:obj.coord_limits(2);
        end;   
        function x = x_ep(obj)
            x = obj.coord_limits(1):obj.coord_stps:(obj.coord_limits(2)+obj.coord_stps);
        end;         
        
    end
end
