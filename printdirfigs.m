function printdirfigs(varargin)


if nargin<1 || isempty(varargin{1})
    files=dir('*.fig');
else
    files=varargin{1};
    validateattributes(files,{'struct','cell'});
end

if nargin<2 || isempty(varargin{2})
    filetype='pdf';
else
    filetype=varargin{2};
end


for f=1:length(files)
    if iscell(files)
        pdfname=files{f};
    else
        pdfname=strsplit(files(f).name,'.');
        pdfname=pdfname{1};
    end
    
    hf = openfig(pdfname,'reuse');
    pos = hf.Position;

    printfig(hf,pdfname,filetype,[],0);
    
    close(hf);
end