function [Basis,Elem] = setBsplineGrid_func(knotvectorU, knotvectorV, knotvectorW,uknotvectorU,uknotvectorV,uknotvectorW, param, maxlevel)

nelemU = param.nelemU;
nelemV = param.nelemV;
nelemW = param.nelemW;
nobU = param.nobU;
nobV = param.nobV;
nobW = param.nobW;
Elem = cell(maxlevel,1);
Basis = cell(maxlevel,1);
sizem = 0;

for level = 1:maxlevel,
    level
    
    knot_cu = knotvectorU{level,1};
    knot_cv = knotvectorV{level,1};
    knot_cw = knotvectorW{level,1};
    
    if(level<=maxlevel-1)
        knot_fu = knotvectorU{level+1,1};
        knot_fv = knotvectorV{level+1,1};
        knot_fw = knotvectorW{level+1,1};
    else
        knot_fu = 0;
        knot_fv = 0;
        knot_fw = 0;
    end
    
    if(level>1)
        Elem{level,7} = parent_elem;
    else
        Elem{level,7} = -1;
    end
    
    [knot_ind, chdElem, Parcell, IEN, cell_centre, supp_i, connectNodes] = storeElemArray_node(level,param,knot_cu, knot_cv, knot_cw, sizem,maxlevel);
    
    parent_elem = Parcell;
    
    %knot indices of each element
    Elem{level,1} = knot_ind;
    
    %activity flag of each element: passive or active
    if(level==1)
        Elem{level,2} = ones(nelemU(level,1)*nelemV(level,1)*nelemW(level,1),1);
    else
        Elem{level,2} = zeros(nelemU(level,1)*nelemV(level,1)*nelemW(level,1),1);
    end
    
    %the non zero splines at the current refinement level over each element
    Elem{level,3} = IEN;
    
    %store nodes associated with each element
    Elem{level,6} = -1;
    
    
    %active cell index
    Elem{level,8} = 0;
    
    %nodes
    Elem{level,9} = connectNodes;
    %store children cells
    if(level<=maxlevel-1)
        Elem{level,4} = chdElem;
    else
        Elem{level,4} = -1*ones(1,8);
    end
    
    %store centroids of each element
    Elem{level,5} = cell_centre;
    
    % Basis Data Structure
    [BB,Chd_Basis,Coeff_Basis,bdyflag] = storeBasisArray(level,param,knot_cu, knot_cv, knot_cw,knot_fu, knot_fv, knot_fw,maxlevel);
    
    %i,j,k index of each B-spline
    Basis{level,1} = BB;
    
    %active flag of the B-spline
    if(level==1)
        Basis{level,2} = ones(nobU(level)*nobV(level)*nobW(level),1);
    else
        Basis{level,2} = zeros(nobU(level)*nobV(level)*nobW(level),1);
    end
    
    %children B-splines
    Basis{level,3} = Chd_Basis;
    
    %coeff matrix of each B-spline
    Basis{level,4} = Coeff_Basis;
    
    %for each B-spline the cell that have the support
    Basis{level,5} = supp_i;
    
    %truncation flag
    Basis{level,6} = zeros(nobU(level)*nobV(level)*nobW(level),1);
    
    %refinement flag
    Basis{level,7} = zeros(nobU(level)*nobV(level)*nobW(level),1);
    %boundary flag
    Basis{level,8} = bdyflag;
    %active B-spline index
    Basis{level,9} = 0;

    sizem = sizem + size(uknotvectorV{level,1},2)*size(uknotvectorU{level,1},2)*size(uknotvectorW{level,1},2);
    %     matFileName = sprintf('mat%d.mat', level);
    %     Bas_t{1,1:end} = Basis{level,1:end};
    %     save(matFileName,'Bas_t');
    %     clear('Basis{level,:}');
end
fieldElem = {'knot_ind','flag_active','IEN','chdElem','cell_centre','node','parElem','actE','nodes'};
Elem = cell2struct(Elem,fieldElem,2);
fieldBasis = {'basis_ind','flag_active','chdBasis','coeffBasis','suppCell','flag_trunc','flag_ref','flag_bdy','actB'};
Basis = cell2struct(Basis,fieldBasis,2);
disp('finish!');
end