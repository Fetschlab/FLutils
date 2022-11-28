function subjIDs = generateSubjectCodes(nSubj,nlet)

% SJ 08-2021
% autogenerate 3-letter codes for human subjects
if nargin<2, nlet = 3; end
if nargin<1, nSubj = 50; end

% allow any letters for all three positions?
alph= 'A':'Z';

for s=1:nSubj
    subjIDs{s,1} = alph(randi(numel(alph),1,nlet));
end



% encryption code?
%om = 'SJJ';
%key = 4;
%ec = char(mod(om + key - 'A', 26) + 'A')
