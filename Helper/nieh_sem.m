function sem1 = nieh_sem(varargin)

if nargin==1
    mat1 = varargin{1};
    sem1 = std(mat1)./sqrt(length(mat1));
    
elseif nargin==2
    
    mat1 = varargin{1};
    dim1 = varargin{2};
    sem1 = std(mat1,0,dim1)./sqrt(size(mat1,dim1));
    
end

