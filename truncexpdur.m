function dur = truncexpdur(meandur,mindur,maxdur,nTrials)

dur = round(exprnd(meandur,nTrials(1),1));

while any(dur<=mindur | dur>=maxdur)
    inds = dur<=mindur | dur>=maxdur;
    dur(inds) = round(exprnd(meandur,sum(inds),1));
end
