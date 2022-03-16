function varargout = load_stl(file)
%
% code borrowed from 'stlread.m'
%

if ~exist(file,'file')
    error(['File ''%s'' not found. If the file is not on MATLAB''s path' ...
        ', be sure to specify the full path to the file.'], file);
end

fid = fopen(file,'r');
if ~isempty(ferror(fid))
    error(lasterror); %#ok
end

M = fread(fid,inf,'uint8=>uint8');
fclose(fid);

if isequal(convertCharsToStrings(char(M(1:5))),'solid')
    [vertices,faces,normals,name] = stlRead(file);
    f = faces;
    v = vertices;
    n = size(f,1);
else
    [f,v,n] = stlbinary(M);
end
%if( isbinary(M) ) % This may not be a reliable test
%    [f,v,n] = stlbinary(M);
%else
%    [f,v,n] = stlascii(M);
%end

varargout = cell(1,nargout);
switch nargout
    case 2
        varargout{1} = f;
        varargout{2} = v;
    case 3
        varargout{1} = f;
        varargout{2} = v;
        varargout{3} = n;
    otherwise
        varargout{1} = struct('faces',f,'vertices',v);
end
end


function [F,V,N] = stlbinary(M)

F = [];
V = [];
N = [];

if length(M) < 84
    error('MATLAB:stlread:incorrectFormat', ...
        'Incomplete header information in binary STL file.');
end

% Bytes 81-84 are an unsigned 32-bit integer specifying the number of faces
% that follow.
numFaces = typecast(M(81:84),'uint32');
if numFaces == 0
    warning('MATLAB:stlread:nodata','No data in STL file.');
    return
end

T = M(85:end);

% numFaces = floor(size(T,1)/50);

F = NaN(numFaces,3);
V = NaN(3*numFaces,3);
N = NaN(numFaces,3);

numRead = 0;
while numRead < numFaces
    % Each facet is 50 bytes
    %  - Three single precision values specifying the face normal vector
    %  - Three single precision values specifying the first vertex (XYZ)
    %  - Three single precision values specifying the second vertex (XYZ)
    %  - Three single precision values specifying the third vertex (XYZ)
    %  - Two unused bytes
    i1    = 50 * numRead + 1;
    i2    = i1 + 50 - 1;
    facet = T(i1:i2)';
    
    n  = typecast(facet(1:12),'single');
    v1 = typecast(facet(13:24),'single');
    v2 = typecast(facet(25:36),'single');
    v3 = typecast(facet(37:48),'single');
    
    n = double(n);
    v = double([v1; v2; v3]);
    
    % Figure out where to fit these new vertices, and the face, in the
    % larger F and V collections.
    fInd  = numRead + 1;
    vInd1 = 3 * (fInd - 1) + 1;
    vInd2 = vInd1 + 3 - 1;
    
    V(vInd1:vInd2,:) = v;
    F(fInd,:)        = vInd1:vInd2;
    N(fInd,:)        = n;
    
    numRead = numRead + 1;
end
end


function [F,V,N] = stlascii(M)
warning('MATLAB:stlread:ascii','ASCII STL files currently not supported.');
F = [];
V = [];
N = [];
end

% TODO: Change the testing criteria! Some binary STL files still begin with
% 'solid'.
function tf = isbinary(A)
% ISBINARY uses the first line of an STL file to identify its format.
if isempty(A) || length(A) < 5
    error('MATLAB:stlread:incorrectFormat', ...
        'File does not appear to be an ASCII or binary STL file.');
end
if strcmpi('solid',char(A(1:5)'))
    tf = false; % ASCII
else
    tf = true;  % Binary
end
end





function [v, f, n, name] = stlRead(fileName)
%STLREAD reads any STL file not depending on its format
%V are the vertices
%F are the faces
%N are the normals
%NAME is the name of the STL object (NOT the name of the STL file)

format = stlGetFormat(fileName);
if strcmp(format,'ascii')
  [v,f,n,name] = stlReadAscii(fileName);
elseif strcmp(format,'binary')
  [v,f,n,name] = stlReadBinary(fileName);
end
end

function format = stlGetFormat(fileName)
%STLGETFORMAT identifies the format of the STL file and returns 'binary' or
%'ascii'

fid = fopen(fileName);
% Check the file size first, since binary files MUST have a size of 84+(50*n)
fseek(fid,0,1);         % Go to the end of the file
fidSIZE = ftell(fid);   % Check the size of the file
if rem(fidSIZE-84,50) > 0
    format = 'ascii';
else
    % Files with a size of 84+(50*n), might be either ascii or binary...
    % Read first 80 characters of the file.
    % For an ASCII file, the data should begin immediately (give or take a few
    % blank lines or spaces) and the first word must be 'solid'.
    % For a binary file, the first 80 characters contains the header.
    % It is bad practice to begin the header of a binary file with the word
    % 'solid', so it can be used to identify whether the file is ASCII or
    % binary.
    fseek(fid,0,-1); % go to the beginning of the file
    header = strtrim(char(fread(fid,80,'uchar')')); % trim leading and trailing spaces
    isSolid = strcmp(header(1:min(5,length(header))),'solid'); % take first 5 char
    fseek(fid,-80,1); % go to the end of the file minus 80 characters
    tail = char(fread(fid,80,'uchar')');
    isEndSolid = findstr(tail,'endsolid');
    
    % Double check by reading the last 80 characters of the file.
    % For an ASCII file, the data should end (give or take a few
    % blank lines or spaces) with 'endsolid <object_name>'.
    % If the last 80 characters contains the word 'endsolid' then this
    % confirms that the file is indeed ASCII.
    if isSolid & isEndSolid
        format = 'ascii';
    else
        format = 'binary';
    end
end
fclose(fid);
end

function [v, f, n, name] = stlReadBinary(fileName)
%STLREADBINARY reads a STL file written in BINARY format
%V are the vertices
%F are the faces
%N are the normals
%NAME is the name of the STL object (NOT the name of the STL file)

%=======================
% STL binary file format
%=======================
% Binary STL files have an 84 byte header followed by 50-byte records, each
% describing a single facet of the mesh.  Technically each facet could be
% any 2D shape, but that would screw up the 50-byte-per-facet structure, so
% in practice only triangular facets are used.  The present code ONLY works
% for meshes composed of triangular facets.
%
% HEADER:
% 80 bytes:  Header text
% 4 bytes:   (int) The number of facets in the STL mesh
%
% DATA:
% 4 bytes:  (float) normal x
% 4 bytes:  (float) normal y
% 4 bytes:  (float) normal z
% 4 bytes:  (float) vertex1 x
% 4 bytes:  (float) vertex1 y
% 4 bytes:  (float) vertex1 z
% 4 bytes:  (float) vertex2 x
% 4 bytes:  (float) vertex2 y
% 4 bytes:  (float) vertex2 z
% 4 bytes:  (float) vertex3 x
% 4 bytes:  (float) vertex3 y
% 4 bytes:  (float) vertex3 z
% 2 bytes:  Padding to make the data for each facet 50-bytes in length
%   ...and repeat for next facet... 

fid = fopen(fileName);
header = fread(fid,80,'int8'); % reading header's 80 bytes
name = deblank(native2unicode(header,'ascii')');
if isempty(name)
    name = 'Unnamed Object'; % no object name in binary files!
end
nfaces = fread(fid,1,'int32');  % reading the number of facets in the stl file (next 4 byters)
nvert = 3*nfaces; % number of vertices
% reserve memory for vectors (increase the processing speed)
n = zeros(nfaces,3);
v = zeros(nvert,3);
f = zeros(nfaces,3);
for i = 1 : nfaces % read the data for each facet
    tmp = fread(fid,3*4,'float'); % read coordinates
    n(i,:) = tmp(1:3); % x,y,z components of the facet's normal vector
    v(3*i-2,:) = tmp(4:6); % x,y,z coordinates of vertex 1
    v(3*i-1,:) = tmp(7:9); % x,y,z coordinates of vertex 2
    v(3*i,:) = tmp(10:12); % x,y,z coordinates of vertex 3
    f(i,:) = [3*i-2 3*i-1 3*i]; % face
    fread(fid,1,'int16'); % Move to the start of the next facet (2 bytes of padding)
end
fclose(fid);
% slim the file (delete duplicated vertices)
[v,f] = stlSlimVerts(v,f);
end


function [vnew, fnew]= stlSlimVerts(v, f)
% PATCHSLIM removes duplicate vertices in surface meshes.
% 
% This function finds and removes duplicate vertices.
%
% USAGE: [v, f]=patchslim(v, f)
%
% Where v is the vertex list and f is the face list specifying vertex
% connectivity.
%
% v contains the vertices for all triangles [3*n x 3].
% f contains the vertex lists defining each triangle face [n x 3].
%
% This will reduce the size of typical v matrix by about a factor of 6.
%
% For more information see:
%  http://www.esmonde-white.com/home/diversions/matlab-program-for-loading-stl-files
%
% Francis Esmonde-White, May 2010

if ~exist('v','var')
    error('The vertex list (v) must be specified.');
end
if ~exist('f','var')
    error('The vertex connectivity of the triangle faces (f) must be specified.');
end

[vnew, indexm, indexn] =  unique(v, 'rows');
fnew = indexn(f);
end


function [v, f, n, name] = stlReadAscii(fileName)
%STLREADASCII reads a STL file written in ASCII format
%V are the vertices
%F are the faces
%N are the normals
%NAME is the name of the STL object (NOT the name of the STL file)

%======================
% STL ascii file format
%======================
% ASCII STL files have the following structure.  Technically each facet
% could be any 2D shape, but in practice only triangular facets tend to be
% used.  The present code ONLY works for meshes composed of triangular
% facets.
%
% solid object_name
% facet normal x y z
%   outer loop
%     vertex x y z
%     vertex x y z
%     vertex x y z
%   endloop
% endfacet
%
% <Repeat for all facets...>
%
% endsolid object_name

fid = fopen(fileName);
cellcontent = textscan(fid,'%s','delimiter','\n'); % read all the file and put content in cells
content = cellcontent{:}(logical(~strcmp(cellcontent{:},''))); % remove all blank lines
fclose(fid);

% read the STL name
line1 = char(content(1));
if (size(line1,2) >= 7)
    name = line1(7:end);
else
    name = 'Unnamed Object';
end

% read the vector normals
normals = char(content(logical(strncmp(content,'facet normal',12))));
n = str2num(normals(:,13:end));

% read the vertex coordinates (vertices)
vertices = char(content(logical(strncmp(content,'vertex',6))));
v = str2num(vertices(:,7:end));
nvert = size(vertices,1); % number of vertices
nfaces = sum(strcmp(content,'endfacet')); % number of faces
if (nvert == 3*nfaces)
    f = reshape(1:nvert,[3 nfaces])'; % create faces
end

% slim the file (delete duplicated vertices)
[v,f] = stlSlimVerts(v,f);
end
