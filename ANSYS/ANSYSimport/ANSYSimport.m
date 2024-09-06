function model = ANSYSimport(filename)
% Imports ANSYS mesh and results data (structural elements and corner nodes only)
%  by M.M. Pedersen, mmp@eng.au.dk
%  Aarhus University, Denmark, 2018
%
% Supported element types:
%   SOLID187: 3D tetrahedral solid
%   SOLID186: 3D hexahedral solid
%   PLANE42,  PLANE82,  PLANE182, PLANE183: 2D 4/8-node quad/triangle
%   SHELL181, SHELL281, SHELL63:            3D 4/8-node shell
%
% INPUTS:
%   [path\]filename for input file (string).
% 
% OUTPUTS:
% model.filename = path + filename of the model file
% model.raw      = raw tables from ANSYS
%   etlist       = element type list: [type_no, ansys_name, internal_type]
%   nlist        = node list: [node_no, x, y, z]
%   elist        = element list: [elem_no, type_no, mat_no, node1, ..., node8, com_x, com_y, com_z]
%   flist        = face list: [node1, node2, node3, node4, elem_no, centr_x, centr_y, centr_z]
%   stress       = stress results: [node_no, Sx, Sy, Sz, Sxy, Syz, Sxz]
%   disp         = displacement results: [node_no, Ux, Uy, Uz, Usum]
% model.surf     = same as above, but containing surface nodes/elements only
%   nlist        = same as above, but containing surface nodes/elements only
%   elist        = same as above, but containing surface nodes/elements only
%   flist        = same as above, but containing surface nodes/elements only
%   stress       = same as above, but containing surface nodes/elements only
%   disp         = same as above, but containing surface nodes/elements only
%   mapping      = maps between names and indices for nodes, elements and faces
%     node2ni    = node_no (externally referenced number) -> node index in nlist
%     elem2ei    = elem_no -> element index in elist
%     face2elem  = face index -> element name
%     ni2face    = node index -> face index
%     ni2ei      = node index -> element index

    hwb = waitbar(0.1,'Reading model file');

    % read text file containing model definition
    raw = read_ansys_file(filename);
    nlist  = raw.nlist;
    elist  = raw.elist;
    etlist = raw.etlist;
    
    node2ni = create_node_map(nlist);
    
    % generate faces for plotting
    waitbar(0.4,hwb,'Generating face table');
    flist = generate_faces(node2ni,elist,etlist);
    raw.flist = flist;

    % general model data
    model.filename = filename;
    model.raw = raw;
    
    % reduction of solid model to surface model
    waitbar(0.6,hwb,'Removing collapsed faces');
    flist = remove_collapsed(flist);
    [nlist,elist,flist] = remap(nlist,elist,flist);

	flist = remove_dup_faces(flist,nlist,hwb);
    [nlist,elist,flist] = remap(nlist,elist,flist);
    
    mapping = create_number_maps(nlist,elist,flist);
    
    % return data for surface model
    model.surf.nlist = nlist;
    model.surf.elist = elist;
    model.surf.flist = flist;
    model.surf.mapping    = mapping;
    model.surf.etypes     = unique(elist(:,2));
           
    % build tables of stress and displacement results for surface model
    surf_nodes = model.surf.nlist(:,1);
    model.surf.stress = model.raw.stress(mapping.node2ni(surf_nodes),:);
    model.surf.disp   = model.raw.disp(mapping.node2ni(surf_nodes),:);
    
    waitbar(1,hwb,'Done.');
    delete(hwb);
    
end


% Sub functions
function raw = read_ansys_file(filename)
    
    % read text file
    fid = fopen(filename);
    if fid==-1
        error('Can''t find or open file: %s\n',filename); 
    end
    file    = textscan(fid,'%s','delimiter',newline,'whitespace','');
    file    = file{1};
    fclose all;
    
    n_lines = size(file,1);

    % initialize output arrays with lots of extra space (trimmed later)
    nlist   = zeros(n_lines,4);     % nodeno, x,y,z coordinates
    elist   = zeros(n_lines,11);    % elemno, typeno, matno, up to 8 corner nodes
    etlist  = zeros(100,1);         % typeno, ansys element name, dim (1=planar, 2=shell, 3=solid)
    stress  = zeros(n_lines,7);     % nodeno, Sx, Sy, Sz, Sxy, Syz, Sxz
    disp    = zeros(n_lines,5);     % nodeno, Ux, Uy, Uz, Usum
    
    % reset counters
    nn  = 0;
    ne  = 0;
    net = 0;
    
    % scan line by line and extract values
    i = 0;
    while i < n_lines

        i = i+1;
        cur_line = file{i};

        
        % read ETLIST
        if strncmp(cur_line,'ET',2)
            
            while ~strncmp(cur_line,'',1)
                typeno = str2double(cur_line(4:5));
                ansys_name = str2double(cur_line(8:11));

                % internal type no, for face generation
                switch ansys_name
                    case {42,82,182,183,181,281,63} % planar/shell
                        dim = 1;
                    case {187,200} % solid tet
                        dim = 2;
                    case {185,186} % solid hex
                        dim = 3;
                    otherwise % unsupported
                        dim = 99; 
                end

                % record
                net = net + 1;
                etlist(net,1) = typeno;
                etlist(net,2) = ansys_name;
                etlist(net,3) = dim;
                
                % goto next line
                i = i+1;
                cur_line = file{i};
            end
            
            % trim used entries
            etlist(net+1:end,:) = [];
            
        end
        
                
        % read NLIST
        if strncmp(cur_line,'   NODE',7)
            
            i = i+1;
            cur_line = file{i};
            
            while ~strncmp(cur_line,'',1)
                % record line
                nn = nn+1;
                nlist(nn,:) = read_node_line(cur_line);
                
                % goto next line
                i = i+1;
                cur_line = file{i};
            end
            
            % trim used entries
            nlist(nn+1:end,:)=[];
            
        end

        
        % read ELIST
        if strncmp(cur_line,'    ELEM',8)
            
            i = i+2;
            cur_line = file{i};
            
            while ~strncmp(cur_line,'',1) 

                if ~strncmp(cur_line,'        ',8) 
                    % record line
                    ne = ne+1;
                    elist(ne,:) = read_elem_line(cur_line,etlist);
                end

                % goto next line
                i = i+1;
                if i<=length(file)
                    cur_line = file{i};
                else
                    break
                end
                
            end
            
            % trim used entries
            elist(ne+1:end,:)=[];
            
        end
        
        
        % read PRNSOL,S
        if strncmp(cur_line,' PRINT S',8)
            
            % skip warning lines etc.
            while ~strncmp(cur_line,'    NODE',8)
                i = i+1;
                cur_line = file{i};
            end
            
            i = i+1;
            cur_line = file{i};
                
            nn = 0;
            
            while ~strncmp(cur_line,'',1)
                % record line
                nn = nn+1;
                stress(nn,:) = read_node_line(cur_line);
                
                % goto next line
                i = i+1;
                cur_line = file{i};
            end
            
            % trim used entries
            stress(nn+1:end,:)=[];
            
        end

        % read PRNSOL,U
        if strncmp(cur_line,' PRINT U',8)
            
            % skip warning lines etc.
            while ~strncmp(cur_line,'    NODE',8)
                i = i+1;
                cur_line = file{i};
            end
            
            i = i+1;
            cur_line = file{i};
            nn = 0;
            
            while ~strncmp(cur_line,'',1)
                % record line
                nn = nn+1;
                disp(nn,:) = read_node_line(cur_line);
                
                % goto next line
                i = i+1;
                if i<=n_lines
                    cur_line = file{i};
                else
                    break
                end
            end
            
            % trim used entries
            disp(nn+1:end,:)=[];
            
        end
        
    end

    % return all results in struct
    raw.etlist  = etlist;
    raw.nlist   = nlist;
    raw.elist   = elist;
    raw.stress  = stress;
    raw.disp    = disp;
      
end

function nlist = read_node_line(cur_line)
% NODE        X                   Y                   Z
%  1    17.6776695300         17.6776695300         0.00000000000     

    str = strsplit(strtrim(cur_line),' ');
    nlist = cellfun(@(str) sscanf(str,'%g'),str);

end

function elist = read_elem_line(cur_line,etlist)
% ELEM MAT TYP REL ESY SEC    NODES
%  1    2   3   4   5   6     7    8     9     10    11    12    13    14
% 316   4   4   1   0   4     74   444   443    75   370   251   252   371 (hex)
%5917   2   2   2   0   2    245   710   711  4244  (shell/planar)
   
    str = strsplit(strtrim(cur_line),' ');
    elem_line = cellfun(@(str) sscanf(str,'%d'),str);

    et = elem_line(3);

    % get no. corner nodes
    if etlist(et,3) == 3 % hex
        npe = 8;
    else % tet/planar/shell 
        npe = 4;
    end
    
    elist = zeros(1,11);
    elist(1) = elem_line(1); % elem no.
    elist(2) = elem_line(2); % material no.
    elist(3) = elem_line(3); % element type no.
    
    elist(4:4+npe-1) = elem_line(7:7+npe-1); % corner nodes

end
  
function flist = remove_dup_faces(flist,nlist,hwb)
% Remove faces with coincident centroids and returns surface-faces only

    waitbar(0.7,hwb,'Calculating face centroids');
    n_faces = size(flist,1);
    c = calc_centroid(flist,nlist);
    waitbar(0.8,hwb,'Removing internal faces');
    
    % find & remove faces with coincident centroids
    [~,ii,~] = unique(c,'rows','stable'); % ii = index of unique values in c
    repi = setdiff(1:n_faces,ii);   % repi = index of one of the repeated values in c
    repc = c(repi,:);
    [~,idx] = setdiff(c,repc,'rows');
    ufaces = flist(idx,:); % surface faces
    ucentr = c(idx,:);          % associated centroids
    flist = [ufaces ucentr];
    
end

function [nlist,elist,flist] = remap(nlist,elist,flist)
% Reduce node, element and face tables so they only include used entities.
% Furthermore remaps nodes to consecutive numbers starting from 1.
% [~,LocB]=ismember(A,B) Locb contains the index in B for each value in A

    face_verts = flist(:,1:4);
    used_node_idx = unique(face_verts);
    %used_elem_idx = unique(flist(:,5));
    used_elem_names = unique(flist(:,5));
    [~,used_elem_idx] = ismember(used_elem_names,elist(:,1));
    
    % remove unused nodes & elements
    nlist = nlist(used_node_idx,:);
    elist = elist(used_elem_idx,:);
    
    % remap vertex entries in face table to point to new node indices
    [~,remapped_verts] = ismember(face_verts,used_node_idx);
    flist = [remapped_verts flist(:,5:end)];
    
end

function node2ni = create_node_map(nlist)
    
    % node name -> index
    nodes = nlist(:,1);
    node2ni  = zeros(max(nodes),1);
    for ni = 1:size(nodes,1) % ni = internal node number
        n = nodes(ni);       % n  = external node name
        node2ni(n) = ni;
    end
    
end

function mapping = create_number_maps(nlist,elist,flist)
% create maps between different number systems

    n_nodes = size(nlist,1);
    n_elems = size(elist,1);
    n_faces = size(flist,1);
    

    % node name -> index (for surface nodes only)
    node2ni = create_node_map(nlist);
    
    % element name -> index
    elems = elist(:,1);
    elem2ei = zeros(max(elems),1,'uint32');
    for ei = 1:size(elems,1) % ei = internal element number
        e = elems(ei);       % e  = external element name
        elem2ei(e) = ei;
    end
    
    
    % face index -> element name
    face2elem = flist(:,5);
    
   
    % node index -> face no.
    faces = flist(:,1:4);
    ni2face = zeros(n_nodes,10,'uint32');% node2face(node) = face
    nf = zeros(n_nodes,1,'uint8');
    
    for f = 1:n_faces
        
        if faces(f,3) == faces(f,4)
            fnodes = 3;
        else
            fnodes = 4;
        end
        
        for i = 1:fnodes
            ni = faces(f,i); % current node
            nf(ni) = nf(ni) + 1; % counter
            ni2face(ni,nf(ni)) = f; % current face
        end
        
    end
    
    
    % map node-idx -> connected elem-idx
    ni2ei = zeros(n_nodes,15,'uint32');
    ne = zeros(n_nodes,1,'uint8');
    
    for ei = 1:n_elems
        for i = 4:11 % node numbers in elist
            n = elist(ei,i); % current node
            if n>0 && n<=length(node2ni)
                ni = node2ni(n);
                if ni>0
                    ne(ni) = ne(ni) + 1; % counter
                    ni2ei(ni,ne(ni)) = ei; % current elem
                end
            end
        end
    end
    
    % return values
    mapping.node2ni   = node2ni;
    mapping.elem2ei   = elem2ei;
    mapping.face2elem = face2elem;
    mapping.ni2face   = ni2face;
    mapping.ni2ei     = ni2ei;
    
end

function flist = generate_faces(node2idx,elist,etlist)
% Generates table of faces for all elements
% incl. map of facenumber to element index

    % mapping of nodes on tet/hex elements to faces
    tetmap  = [1 2 3 3; 1 2 4 4; 1 3 4 4; 2 3 4 4]; % use 4 vertices for tet faces also.
    hexmap  = [1 2 3 4; 1 2 6 5; 2 3 7 6; 3 7 8 4; 1 4 8 5; 5 6 7 8];
    quadmap = [1 2 3 4];

    % face table: [node_idx1, node_idx2, node_idx3, node_idx4, elem_name]
    n_elems = size(elist,1);
    flist = zeros(6*n_elems,5);
    f = 0;
    
    for ei = 1:n_elems

        % get element type specifics
        cur_elem_type = elist(ei,3);
        elem_name     = elist(ei,1);
        
        switch etlist(cur_elem_type,3)
            case 1 % planar/shell
                nodes   = elist(ei,4:7);
                fmap    = quadmap;
                n_faces = 1;
            case 2 % tet
                nodes   = elist(ei,4:7);
                fmap    = tetmap;
                n_faces = 4;
            case 3 % hex
                nodes   = elist(ei,4:11);
                fmap    = hexmap;
                n_faces = 6;
            otherwise
                n_faces = 0;
        end
        
        % build face table
        for fi = 1:n_faces
            f = f + 1;
            flist(f,1:4) = node2idx(nodes(fmap(fi,:)))';
            flist(f,5)   = elem_name;
        end
        
    end
    
    % remove unused entries
    flist(f+1:end,:) = [];
    
end

function flist = remove_collapsed(flist)
% Remove collapsed faces, i.e. faces where multiple vertex numbers 
% refer to the same nodes. E.g. face = point or face = line.

    n_faces = size(flist,1);
    unique_nodes_on_face = zeros(n_faces,1,'uint32');
    
    for f = 1:n_faces
        
        unique_nodes_on_face(f) = length(unique(flist(f,1:4)));
        
    end
    
    collapsed = [find(unique_nodes_on_face==1); 
    find(unique_nodes_on_face==2)];
    flist(collapsed,:) = [];
    
end

function c = calc_centroid(flist,nlist)

    n_faces = size(flist,1);
    c   = zeros(n_faces,3,'single');
    
    % calculate centroid for all faces
    for f = 1:n_faces
        
        cur_face = unique(flist(f,1:4));
        n_nodes = length(cur_face);
        
        for no = 1:n_nodes
            p = nlist(cur_face(no),2:4);
            c(f,:) = c(f,:) + p;
        end
        
        c(f,:) = c(f,:)/n_nodes;
    end
    
end
